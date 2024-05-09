from pytket.circuit import Circuit, OpType, CircBox, PhasePolyBox, Conditional
from pytket.passes import CustomPass, RepeatWithMetricPass
from pytket.predicates import GateSetPredicate
from pytket.unit_id import Bit, Qubit


fswap_circ = Circuit(2, name="FSWAP").CZ(0, 1).SWAP(0, 1)

fswap = CircBox(fswap_circ)


def gadgetise_hadamards(circ: Circuit) -> Circuit:
    h_count = circ.n_gates_of_type(OpType.H)

    circ_prime = Circuit(circ.n_qubits)
    z_ancillas = circ_prime.add_q_register("z_ancillas", h_count)
    ancilla_bits = circ_prime.add_c_register("bits", h_count)

    for ancilla in z_ancillas:
        circ_prime.H(ancilla)

    circ_prime.add_barrier(list(z_ancillas))

    ancilla_index = 0

    for cmd in circ:
        if cmd.op.type != OpType.H:
            circ_prime.add_gate(cmd.op, cmd.args)
        else:
            circ_prime.add_gate(fswap, [cmd.qubits[0], z_ancillas[ancilla_index]])
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


def reverse_circuit(circ: Circuit) -> Circuit:
    new_circ = _initialise_registers(circ)

    for cmd in reversed(circ.get_commands()):
        if cmd.op.type == OpType.Barrier:
            new_circ.add_barrier(cmd.qubits)
        else:
            new_circ.add_gate(cmd.op, cmd.args)

    return new_circ


pauli_prop_predicate = GateSetPredicate(
    {OpType.Measure, OpType.CircBox, OpType.PhasePolyBox, OpType.Conditional}
)


def _get_terminal_pauli_x_args(circ: Circuit) -> tuple[Bit, Qubit]:
    for cmd in reversed(circ.get_commands()):
        if cmd.op.type == OpType.Conditional:
            if cmd.op.op.type == OpType.X:
                x_pauli_bit = cmd.args[0]
                x_pauli_qubit = cmd.args[1]
    return x_pauli_bit, x_pauli_qubit


def phasepolybox_to_conjugation(poly_box: PhasePolyBox, x_index: int) -> CircBox:
    poly_circ = poly_box.get_circuit()
    poly_circ_dg = poly_circ.dagger()
    poly_circ.X(x_index)
    poly_circ.append(poly_circ_dg)
    poly_circ.name = "U X U†"
    return CircBox(poly_circ)


def propogate_terminal_pauli_x_gate(circ: Circuit) -> Circuit:
    reversed_circ = reverse_circuit(circ)
    pauli_prop_predicate.verify(reversed_circ)
    circ_prime = _initialise_registers(reversed_circ)
    pauli_x_args = _get_terminal_pauli_x_args(reversed_circ)
    found_match = False
    for cmd in reversed_circ:
        if cmd.op.type == OpType.PhasePolyBox:
            if pauli_x_args[1] in cmd.qubits and not found_match:
                found_match = True
                uxudg_box = phasepolybox_to_conjugation(
                    cmd.op, pauli_x_args[1].index[0]
                )
                circ_prime.add_gate(
                    uxudg_box,
                    cmd.qubits,
                    condition_bits=[pauli_x_args[0]],
                    condition_value=1,
                )
            else:
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
            # elif cmd.op.circuit_name == "FSWAP":
            #    pass
            else:
                circ_prime.add_gate(cmd.op, cmd.args)

        elif cmd.op.type == OpType.Barrier:
            circ_prime.add_barrier(cmd.qubits)
        else:
            circ_prime.add_gate(cmd.op.type, cmd.args)

    return reverse_circuit(circ_prime)


propogate_terminal_pauli = CustomPass(propogate_terminal_pauli_x_gate)


def is_conditional_pauli_x(operation: Conditional) -> bool:
    return operation.op.type == OpType.X


def get_n_conditional_paulis(circ: Circuit) -> int:
    conditional_ops = circ.ops_of_type(OpType.Conditional)
    conditional_xs = list(filter(is_conditional_pauli_x, conditional_ops))
    return len(conditional_xs)


propogate_all_terminal_paulis = RepeatWithMetricPass(
    propogate_terminal_pauli, get_n_conditional_paulis
)
