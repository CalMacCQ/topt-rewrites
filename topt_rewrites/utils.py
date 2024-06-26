from pytket._tket.circuit import Circuit, Conditional, OpType, PhasePolyBox
from pytket.predicates import NoSymbolsPredicate


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
    """Check that the underlying Circuit for a PhasePolyBox is Clifford + T."""
    circ = ppb.get_circuit()
    return check_rz_angles(circ)


def _is_conditional_pauli_x(operation: Conditional) -> bool:
    return operation.op.type == OpType.X


def get_n_conditional_paulis(circ: Circuit) -> int:
    """Return the number of Conditonal-X gates in a Circuit."""
    conditional_ops = circ.ops_of_type(OpType.Conditional)
    conditional_xs = list(filter(_is_conditional_pauli_x, conditional_ops))
    return len(conditional_xs)


def initialise_registers(circ: Circuit) -> Circuit:
    circ_prime = Circuit()
    for qb in circ.qubits:
        circ_prime.add_qubit(qb)

    for bit in circ.bits:
        circ_prime.add_bit(bit)
    return circ_prime


def reverse_circuit(circ: Circuit) -> Circuit:
    new_circ = initialise_registers(circ)

    for cmd in reversed(circ.get_commands()):
        if cmd.op.type == OpType.Barrier:
            new_circ.add_barrier(cmd.qubits)
        else:
            new_circ.add_gate(cmd.op, cmd.args)

    return new_circ
