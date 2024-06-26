"""Circuit preprocessing for T gate optimisation."""  # noqa: INP001

from __future__ import annotations

from pytket.circuit import (
    CircBox,
    Circuit,
    Conditional,
    OpType,
    PhasePolyBox,
)

from pytket.predicates import GateSetPredicate
from pytket.unit_id import Bit, Qubit  # noqa: TCH002

from .utils import (
    initialise_registers,
    reverse_circuit,
    get_n_conditional_paulis,
)


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


# PROPOGATE_ALL_TERMINAL_PAULIS = RepeatWithMetricPass(
#    PROPAGATE_TERMINAL_PAULI,
#    get_n_conditional_paulis,
# )


def _get_v_box(circ: Circuit) -> CircBox:

    v_circ = initialise_registers(circ)
    v_circ.name = "V"

    reversed_circ = reverse_circuit(circ)

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

    new_circ = initialise_registers(circ)

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
