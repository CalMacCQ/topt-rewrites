"""Circuit preprocessing for T gate optimisation."""  # noqa: INP001

from __future__ import annotations

from pytket.circuit import (
    CircBox,
    Circuit,
    Command,
    Conditional,
    OpType,
    PauliExpCommutingSetBox,
    PhasePolyBox,
)
from pytket.extensions.qiskit import qiskit_to_tk
from pytket.passes import (
    CustomPass,
    DecomposeBoxes,
    RepeatWithMetricPass,
)
from pytket.pauli import Pauli, QubitPauliTensor
from pytket.predicates import GateSetPredicate, NoSymbolsPredicate
from pytket.tableau import UnitaryTableau, UnitaryTableauBox
from pytket.unit_id import Bit, Qubit  # noqa: TCH002
from qiskit.synthesis import synth_cnot_count_full_pmh

FSWAP_CIRC = Circuit(2, name="FSWAP").CZ(0, 1).SWAP(0, 1)

FSWAP = CircBox(FSWAP_CIRC)

HADAMARD_REPLACE_PREDICATE = GateSetPredicate({OpType.H, OpType.PhasePolyBox})


def get_n_internal_hadamards(circ: Circuit) -> int:
    """Return the number of Hadamard gates between the first and last non-Clifford gate in the Circuit."""
    if not HADAMARD_REPLACE_PREDICATE.verify(circ):
        pred_msg = "Circuit must contain only OpType.H and OpType.PhasePolyBox OpTypes."
        raise ValueError(pred_msg)

    total_h_count = circ.n_gates_of_type(OpType.H)

    # Count number of Hadamards until we encounter a PhasePolyBox
    lhs_count = _count_hadamards(circ.get_commands())

    # Same but from the end of the circuit
    rhs_count = _count_hadamards(reversed(circ.get_commands()))

    # TODO add handling for when a PhasePolyBox is Clifford

    return total_h_count - (lhs_count + rhs_count)


def _count_hadamards(commands: list[Command]) -> int:
    h_count = 0
    for cmd in commands:
        if cmd.op.type == OpType.H:
            h_count += 1
        elif cmd.op.type == OpType.PhasePolyBox:
            break
    return h_count


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


# TODO Refactor this awful propagate_terminal_pauli_x_gate function.

# def propagate_terminal_pauli_x_gate(circ: Circuit) -> Circuit:  # noqa: PLR0912
#    """Propogates a single Pauli X gate to the end of the Circuit."""
#    if not PAULI_PROP_PREDICATE.verify(circ):
#        msg = f"Circuit must be in the {PAULI_PROP_GATES} gateset."
#        raise ValueError(msg)
#    reversed_circ = _reverse_circuit(circ)
#    circ_prime = _initialise_registers(reversed_circ)
#    pauli_x_args = _get_terminal_pauli_x_args(reversed_circ)
#    found_match = False
#    for cmd in reversed_circ:
#        if cmd.op.type == OpType.PhasePolyBox:
#            if pauli_x_args[1] in cmd.qubits and not found_match:
#                found_match = True
#                uxudg_box = _get_conjugation(
#                    cmd.op,
#                    pauli_x_args[1].index[0],
#                )
#                circ_prime.add_gate(
#                    uxudg_box,
#                    cmd.qubits,
#                    condition_bits=[pauli_x_args[0]],
#                    condition_value=1,
#                )
#            # Add PhasePolyBox as usual in both cases
#            circ_prime.add_gate(cmd.op, cmd.qubits)
#
#        elif cmd.op.type == OpType.Measure:
#            circ_prime.Measure(cmd.args[0], cmd.args[1])
#
#        elif cmd.op.type == OpType.Conditional:
#            if cmd.op.op.type != OpType.X:
#                circ_prime.add_gate(cmd.op, cmd.args)
#            elif cmd.op.op.type == OpType.X:
#                if (
#                    tuple(cmd.args) == pauli_x_args
#                ):  # check for matching qubits and bits
#                    pass
#                else:
#                    circ_prime.add_gate(cmd.op, cmd.args)
#            else:
#                circ_prime.add_gate(cmd.op, cmd.args)
#
#        elif cmd.op.type == OpType.CircBox:
#            if cmd.op.circuit_name == "U X U†":
#                pass  # TODO propogate The X through the U X Udg and FSWAP
#            # elif cmd.op.circuit_name == "FSWAP":  # noqa: ERA001
#            #    pass
#            else:
#                circ_prime.add_gate(cmd.op, cmd.args)
#
#        elif cmd.op.type == OpType.Barrier:
#            circ_prime.add_barrier(cmd.qubits)
#        else:
#            circ_prime.add_gate(cmd.op.type, cmd.args)
#
#    return _reverse_circuit(circ_prime)


# PROPAGATE_TERMINAL_PAULI = CustomPass(propagate_terminal_pauli_x_gate)


def _is_conditional_pauli_x(operation: Conditional) -> bool:
    return operation.op.type == OpType.X


def get_n_conditional_paulis(circ: Circuit) -> int:
    """Return the number of Conditonal-X gates in a Circuit."""
    conditional_ops = circ.ops_of_type(OpType.Conditional)
    conditional_xs = list(filter(_is_conditional_pauli_x, conditional_ops))
    return len(conditional_xs)


# PROPOGATE_ALL_TERMINAL_PAULIS = RepeatWithMetricPass(
#    PROPAGATE_TERMINAL_PAULI,
#    get_n_conditional_paulis,
# )


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


def _construct_partial_circuit(circ: Circuit) -> Circuit:
    if get_n_conditional_paulis(circ) < 1:
        msg = "Circuit contains no conditional X gates."
        raise ValueError(msg)

    new_circ = _initialise_registers(circ)

    for cmd in circ:
        if cmd.op.type == Conditional:
            if cmd.op.op.type == OpType.X:
                break
        else:
            new_circ.add_gate(cmd.op, cmd.qubits)
    return new_circ


def _construct_full_circuit(
    circ: Circuit,
    v_operator: CircBox,
    x_index: int,
) -> Circuit:
    normal_circ = _construct_partial_circuit(circ)
    vxvdg_box = _get_conjugation(v_operator, x_index)
    normal_circ.add_gate(
        vxvdg_box,
        circ.qubits,
        condition_bits=[circ.bits[0]],
        condition_value=1,
    )
    return normal_circ


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

PAULI_DICT = {
    Pauli.I: OpType.noop,
    Pauli.X: OpType.X,
    Pauli.Y: OpType.Y,
    Pauli.Z: OpType.Z,
}


def _get_new_pauli(circ: Circuit, tableau: UnitaryTableauBox) -> Circuit:
    pass


def _get_reversible_tableau(pbox: PhasePolyBox) -> UnitaryTableau:
    # cheat by synthesising the CNOT circuit with qiskit and converting
    qc = synth_cnot_count_full_pmh(pbox.linear_transformation, section_size=2)
    qc2 = qc.reverse_bits()  # correct for endianness
    tkc_cnot = qiskit_to_tk(qc2)
    return UnitaryTableau(tkc_cnot)


def _parities_to_pauli_tensors(pbox: PhasePolyBox) -> list[QubitPauliTensor]:
    phase_poly = pbox.phase_polynomial
    tensor_list = []
    for parity, phase in phase_poly.items():
        pauli_list = []
        qubit_list = []
        for count, boolean in enumerate(parity):
            qubit_list.append(Qubit(count))
            if boolean:
                pauli_list.append(Pauli.Z)
            else:
                pauli_list.append(Pauli.I)

        pauli_tensor = QubitPauliTensor(
            qubits=qubit_list,
            paulis=pauli_list,
            coeff=phase,
        )
        tensor_list.append(pauli_tensor)

    return tensor_list


def _get_circuit_fragments(
    pbox: PhasePolyBox,
    input_pauli: QubitPauliTensor,
) -> tuple[Circuit, Circuit]:
    # p_box = S L
    # S: Sequence of {Z, I} Phase gadgets
    # L: Linear reversible transformation

    # Compute Tableau form of L.
    l_tableau: UnitaryTableau = _get_reversible_tableau(pbox)

    # Updated Pauli Tensor P'. P' = L P L†.
    new_pauli: QubitPauliTensor = l_tableau.get_row_product(input_pauli)

    # Create new circuit implementing P'
    pauli_prime_circ = Circuit(pbox.n_qubits)
    for qubit, pauli_op in new_pauli.string.map.items():
        pauli_prime_circ.add_gate(PAULI_DICT[pauli_op], [qubit])

    # Convert (True, False) terms to Pauli.Z, Pauli.I in a QubitPauliTensor.
    pauli_tensors = _parities_to_pauli_tensors(pbox)

    # Get updated phase gadget Sequence S' and synthesise circuit.
    s_prime = _get_updated_paulis(pauli_tensors, new_pauli)
    s_prime_circ = _get_phase_gadget_circuit(s_prime)

    return pauli_prime_circ, s_prime_circ


# def get_clifford_circuit()


def _get_updated_paulis(
    pauli_tensors: list[QubitPauliTensor],
    new_pauli: QubitPauliTensor,
) -> list[QubitPauliTensor]:

    for pauli_op in pauli_tensors:
        if not pauli_op.commutes_with(new_pauli):
            pauli_op.coeff *= 2

    return pauli_tensors


def _get_phase_gadget_circuit(pauli_tensors: list[QubitPauliTensor]) -> Circuit:

    pauli_ops = []
    for tensor in pauli_tensors:
        pauli_list = list(tensor.string.map.values())
        pair = (pauli_list, tensor.coeff.real)
        pauli_ops.append(pair)

    pauli_gadgets_sequence = PauliExpCommutingSetBox(pauli_ops)
    n_qubits = pauli_gadgets_sequence.n_qubits

    pauli_gadget_circ = Circuit(n_qubits).add_gate(
        pauli_gadgets_sequence,
        list(range(n_qubits)),
    )

    DecomposeBoxes().apply(pauli_gadget_circ)

    return pauli_gadget_circ
