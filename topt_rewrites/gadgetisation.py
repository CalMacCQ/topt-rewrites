from __future__ import annotations

from pytket._tket.circuit import CircBox, Circuit, Command, OpType
from pytket.passes import CustomPass
from pytket.predicates import GateSetPredicate

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
