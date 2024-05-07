from pytket.circuit import Circuit, OpType, CircBox
from pytket.passes import CustomPass


new_circ = Circuit(2, name="H-Gadget").H(1).CZ(0, 1).SWAP(0, 1)

H_gadget = CircBox(new_circ)


def gadgetise_hadamards(circ: Circuit) -> Circuit:
    h_count = circ.n_gates_of_type(OpType.H)

    circ_prime = Circuit(circ.n_qubits)
    z_ancillas = circ_prime.add_q_register("z_ancillas", h_count)
    ancilla_bits = circ_prime.add_c_register("bits", h_count)

    ancilla_index = 0

    for cmd in circ:
        if cmd.op.type != OpType.H:
            circ_prime.add_gate(cmd.op.type, cmd.args)
        else:
            circ_prime.add_gate(H_gadget, [cmd.qubits[0], z_ancillas[ancilla_index]])
            circ_prime.Measure(z_ancillas[ancilla_index], ancilla_bits[ancilla_index])
            circ_prime.X(cmd.qubits[0], condition_bits=[ancilla_bits[ancilla_index]], condition_value=1)
            ancilla_index += 1
            
    return circ_prime


replace_hadamards = CustomPass(gadgetise_hadamards)