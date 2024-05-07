from pytket.circuit import Circuit, OpType, CircBox, PhasePolyBox
from pytket.passes import CustomPass
from pytket.predicates import GateSetPredicate
from pytket.unit_id import Bit, Qubit


new_circ = Circuit(2, name="H-Gadget").H(1).CZ(0, 1).SWAP(0, 1)

H_gadget = CircBox(new_circ)


def phasepolybox_to_conjugation(
    poly_box: PhasePolyBox, x_index: int, count: int
) -> CircBox:
    poly_circ = poly_box.get_circuit()
    poly_circ.X(x_index)
    poly_circ.append(poly_circ.dagger())
    poly_circ.name = f"U_{count}X U†"
    return CircBox(poly_circ)


def gadgetise_hadamards(circ: Circuit) -> Circuit:
    h_count = circ.n_gates_of_type(OpType.H)

    circ_prime = Circuit(circ.n_qubits)
    z_ancillas = circ_prime.add_q_register("z_ancillas", h_count)
    ancilla_bits = circ_prime.add_c_register("bits", h_count)

    ancilla_index = 0

    for cmd in circ:
        if cmd.op.type != OpType.H:
            if cmd.op.is_gate():
                circ_prime.add_gate(cmd.op.type, cmd.op.params, cmd.qubits)
            else:
                # Handle CircBox and other boxes
                circ_prime.add_gate(cmd.op, cmd.qubits)
        else:
            circ_prime.add_gate(H_gadget, [cmd.qubits[0], z_ancillas[ancilla_index]])
            circ_prime.Measure(z_ancillas[ancilla_index], ancilla_bits[ancilla_index])
            circ_prime.X(
                cmd.qubits[0],
                condition_bits=[ancilla_bits[ancilla_index]],
                condition_value=1,
            )
            ancilla_index += 1
    return circ_prime


replace_hadamards = CustomPass(gadgetise_hadamards)


def _initialise_registers(circ: Circuit) -> Circuit:
    circ_prime = Circuit()
    for qreg in circ.q_registers:
        for qubit in qreg:
            circ_prime.add_qubit(qubit)

    for creg in circ.c_registers:
        for bit in creg:
            circ_prime.add_bit(bit)
    return circ_prime


pauli_prop_predicate = GateSetPredicate(
    {OpType.Measure, OpType.CircBox, OpType.PhasePolyBox, OpType.Conditional}
)


def get_last_pauli_x_args(circ: Circuit) -> tuple[Qubit, Bit]:
    for cmd in reversed(circ.get_commands()):
        if cmd.op.type == OpType.Conditional:
            x_pauli_qubit = cmd.args[0]
            x_pauli_bit = cmd.args[1]
    return x_pauli_qubit, x_pauli_bit


def propogate_final_pauli_x_gate(circ: Circuit) -> Circuit:
    pauli_prop_predicate.verify(circ)
    circ_prime = _initialise_registers(circ)
    pauli_x_args = get_pauli_x_args(circ)
    phase_poly_count = 0
    for cmd in reversed(circ.get_commands()):
        if cmd.op.type == OpType.Conditional:
            pass
        elif cmd.op.type == PhasePolyBox:
            poly_circ = cmd.op.get_circuit()
            poly_circ.name = "U_f"
            poly_circ.X(0)
            poly_circ.append(poly_circ.dagger())

        elif cmd.op.type == OpType.CircBox:
            circ_prime.add_gate(cmd.op.type, cmd.op.params, cmd.qubits)
        # elif cmd.op.type == OpType.PhasePolyBox:
        #    if
        #        phase_poly_circ = cmd.op.get_circuit()
        #        dg_phase_poly = phase_poly_circ.dagger()
        #    clifford

    return circ_prime
