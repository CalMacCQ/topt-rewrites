"""Circuit preprocessing for T gate optimisation."""  # noqa: INP001

from __future__ import annotations

from pytket.circuit import CircBox, Circuit, Conditional, OpType, PhasePolyBox
from pytket.passes import (
    CustomPass,
    RepeatWithMetricPass,
    ComposePhasePolyBoxes,
    DecomposeBoxes,
    FlattenRegisters,
)
from pytket.predicates import GateSetPredicate, NoSymbolsPredicate
from pytket.unit_id import Bit, Qubit  # noqa: TCH002

FSWAP_CIRC = Circuit(2, name="FSWAP").CZ(0, 1).SWAP(0, 1)

FSWAP = CircBox(FSWAP_CIRC)

HADAMARD_REPLACE_PREDICATE = GateSetPredicate({OpType.H, OpType.PhasePolyBox})


def get_n_internal_hadamards(circ: Circuit) -> int:
    """Return the number of Hadamard gates between the first and last non-Clifford gate in the Circuit."""
    if not HADAMARD_REPLACE_PREDICATE.verify(circ):
        pred_msg = "Circuit must contain only OpType.H and OpType.PhasePolyBox OpTypes."
        raise ValueError(pred_msg)

    total_h_count = circ.n_gates_of_type(OpType.H)

    external_h_count = 0

    for cmd in circ.get_commands():
        if cmd.op.type == OpType.H:
            external_h_count += 1
        elif cmd.op.type == OpType.PhasePolyBox:
            break

    for cmd in reversed(circ.get_commands()):
        if cmd.op.type == OpType.H:
            external_h_count += 1
        elif cmd.op.type == OpType.PhasePolyBox:
            break

    return total_h_count - external_h_count


def gadgetise_hadamards(circ: Circuit) -> Circuit:
    """Replace all Hadamard gates with measurement gadgets."""
    internal_h_count = get_n_internal_hadamards(circ)

    circ_prime = Circuit(circ.n_qubits)
    z_ancillas = circ_prime.add_q_register("z_ancillas", internal_h_count)
    ancilla_bits = circ_prime.add_c_register("bits", internal_h_count)

    for ancilla in z_ancillas:
        circ_prime.H(ancilla)

    circ_prime.add_barrier(list(z_ancillas))

    ancilla_index = 0

    for cmd in circ:
        if cmd.op.type != OpType.H:
            circ_prime.add_gate(cmd.op, cmd.args)
        else:
            circ_prime.add_gate(FSWAP, [cmd.qubits[0], z_ancillas[ancilla_index]])
            circ_prime.Measure(z_ancillas[ancilla_index], ancilla_bits[ancilla_index])
            circ_prime.X(
                cmd.qubits[0],
                condition_bits=[ancilla_bits[ancilla_index]],
                condition_value=1,
            )
            ancilla_index += 1
    return circ_prime


REPLACE_HADAMARDS = CustomPass(gadgetise_hadamards)


def check_rz_angles(circ: Circuit) -> bool:
    """Check that all Rz gates in a Circuit can be implemented with Clifford+T gates."""
    if not NoSymbolsPredicate().verify(circ):
        symbol_msg = "Circuit contains symbolic angles."
        raise ValueError(symbol_msg)

    if circ.n_gates_of_type(OpType.Rz) == 0:
        no_rz_error = "Circuit does not contain any Rz gates."
        raise ValueError(no_rz_error)

    rz_op_list = circ.ops_of_type(OpType.Rz)

    allowed_non_clifford_angles = [0.25, 0.75, 1.25, 1.75]

    for op in rz_op_list:
        if not op.is_clifford():
            if abs(op.params[0]) % 2 in allowed_non_clifford_angles:
                pass
            else:
                return False

    return True


def check_phasepolybox(ppb: PhasePolyBox) -> bool:
    """Check that the underlying Circuit for a PhasePolyBox is Clifford."""
    circ = ppb.get_circuit()
    return check_rz_angles(circ)


def _initialise_registers(circ: Circuit) -> Circuit:
    circ_prime = Circuit()
    for qb in circ.qubits:
        circ_prime.add_qubit(qb)

    for bit in circ.bits:
        circ_prime.add_bit(bit)
    return circ_prime


def _reverse_circuit(circ: Circuit) -> Circuit:
    new_circ = _initialise_registers(circ)

    for cmd in reversed(circ.get_commands()):
        if cmd.op.type == OpType.Barrier:
            new_circ.add_barrier(cmd.qubits)
        else:
            new_circ.add_gate(cmd.op, cmd.args)

    return new_circ


PAULI_PROP_GATES = {
    OpType.PhasePolyBox,
    OpType.H,
    OpType.X,
    OpType.CircBox,
    OpType.Measure,
    OpType.Conditional,
    OpType.Barrier,
}

PAULI_PROP_PREDICATE = GateSetPredicate(PAULI_PROP_GATES)


def _get_terminal_pauli_x_args(circ: Circuit) -> tuple[Bit, Qubit]:
    if get_n_conditional_paulis(circ) < 1:
        msg = "Circuit contains no conditional X gates."
        raise ValueError(msg)
    for cmd in reversed(circ.get_commands()):
        if cmd.op.type == OpType.Conditional and cmd.op.op.type == OpType.X:
            x_pauli_bit = cmd.args[0]
            x_pauli_qubit = cmd.args[1]
    return x_pauli_bit, x_pauli_qubit


def _get_conjugation(box: PhasePolyBox | CircBox, x_index: int) -> CircBox:
    circ = box.get_circuit()
    circ_dg = circ.dagger()
    circ.X(x_index)
    circ.append(circ_dg)
    if isinstance(box, PhasePolyBox):
        circ.name = "U X U†"
    else:
        circ.name = "V X V†"
    return CircBox(circ)


def _get_n_terminal_boxes(circ: Circuit) -> int:
    PAULI_PROP_PREDICATE.verify(circ)
    phase_poly_count = 0
    backwards_circuit = _reverse_circuit(circ)
    for cmd in backwards_circuit:
        match cmd.op.type:
            case OpType.PhasePolyBox:
                phase_poly_count += 1
            case OpType.CircBox if cmd.op.circuit_name == "U X U†":
                phase_poly_count += 1
            case OpType.Conditional if cmd.op.op.type == OpType.X:
                break
    return phase_poly_count


def propagate_terminal_pauli_x_gate(circ: Circuit) -> Circuit:  # noqa: PLR0912
    """Propogates a single Pauli X gate to the end of the Circuit."""
    if not PAULI_PROP_PREDICATE.verify(circ):
        msg = f"Circuit must be in the {PAULI_PROP_GATES} gateset."
        raise ValueError(msg)
    reversed_circ = _reverse_circuit(circ)
    circ_prime = _initialise_registers(reversed_circ)
    pauli_x_args = _get_terminal_pauli_x_args(reversed_circ)
    found_match = False
    for cmd in reversed_circ:
        if cmd.op.type == OpType.PhasePolyBox:
            if pauli_x_args[1] in cmd.qubits and not found_match:
                found_match = True
                uxudg_box = _get_conjugation(
                    cmd.op,
                    pauli_x_args[1].index[0],
                )
                circ_prime.add_gate(
                    uxudg_box,
                    cmd.qubits,
                    condition_bits=[pauli_x_args[0]],
                    condition_value=1,
                )
            # Add PhasePolyBox as usual in both cases
            circ_prime.add_gate(cmd.op, cmd.qubits)

        elif cmd.op.type == OpType.Measure:
            circ_prime.Measure(cmd.args[0], cmd.args[1])

        elif cmd.op.type == OpType.Conditional:
            if cmd.op.op.type != OpType.X:
                circ_prime.add_gate(cmd.op, cmd.args)
            elif cmd.op.op.type == OpType.X:
                if (
                    tuple(cmd.args) == pauli_x_args
                ):  # check for matching qubits and bits
                    pass
                else:
                    circ_prime.add_gate(cmd.op, cmd.args)
            else:
                circ_prime.add_gate(cmd.op, cmd.args)

        elif cmd.op.type == OpType.CircBox:
            if cmd.op.circuit_name == "U X U†":
                pass  # TODO propogate The X through the U X Udg and FSWAP
            # elif cmd.op.circuit_name == "FSWAP":  # noqa: ERA001
            #    pass
            else:
                circ_prime.add_gate(cmd.op, cmd.args)

        elif cmd.op.type == OpType.Barrier:
            circ_prime.add_barrier(cmd.qubits)
        else:
            circ_prime.add_gate(cmd.op.type, cmd.args)

    return _reverse_circuit(circ_prime)


PROPAGATE_TERMINAL_PAULI = CustomPass(propagate_terminal_pauli_x_gate)


def _is_conditional_pauli_x(operation: Conditional) -> bool:
    return operation.op.type == OpType.X


def get_n_conditional_paulis(circ: Circuit) -> int:
    """Return the number of Conditonal-X gates in a Circuit."""
    conditional_ops = circ.ops_of_type(OpType.Conditional)
    conditional_xs = list(filter(_is_conditional_pauli_x, conditional_ops))
    return len(conditional_xs)


PROPOGATE_ALL_TERMINAL_PAULIS = RepeatWithMetricPass(
    PROPAGATE_TERMINAL_PAULI,
    get_n_conditional_paulis,
)


def _get_v_box(circ: Circuit) -> CircBox:

    v_circ = _initialise_registers(circ)
    v_circ.name = "V"

    reversed_circ = _reverse_circuit(circ)

    for cmd in reversed_circ:
        if cmd.op.type == OpType.Measure:
            pass
        elif cmd.op.type == OpType.Conditional:
            if cmd.op.op.type == OpType.X:
                break
        else:
            v_circ.add_gate(cmd.op, cmd.qubits)

    # Remove unnecessary classical bits.
    v_circ.remove_blank_wires()

    return CircBox(v_circ)


# from pytket.circuit.display import view_browser as draw
#
#
# def _build_unitary_subcircuit(circ: Circuit) -> Circuit:
#    circ_prime = _initialise_registers(circ)
#
#    for cmd in circ:
#        if cmd.op.type == OpType.Measure or OpType.Conditional:
#            pass
#        elif cmd.op.type == OpType.Barrier:
#            circ_prime.add_barrier(cmd.qubits)
#        else:
#            circ_prime.add_gate(cmd.op, cmd.args)
#
#        circ_prime.remove_blank_wires()
#        DecomposeBoxes().apply(circ_prime)
#        draw(circ_prime)
#        FlattenRegisters().apply(circ_prime)
#        ComposePhasePolyBoxes().apply(circ_prime)
#        circ_prime.add_gate(cmd.op, cmd.args)
#    return circ_prime
#
