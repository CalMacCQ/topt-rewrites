from pytket.circuit import Circuit, OpType, PhasePolyBox
from pytket.circuit.display import view_browser as draw
from pytket.passes import ComposePhasePolyBoxes, DecomposeBoxes
from pytket.pauli import QubitPauliTensor, Pauli
from pytket.qasm import circuit_to_qasm, circuit_from_qasm

from topt_rewrites.main import (
    # PROPAGATE_TERMINAL_PAULI,
    REPLACE_HADAMARDS,
    check_rz_angles,
    get_n_conditional_paulis,
    get_n_internal_hadamards,
    _get_circuit_fragments,
)

# _get_cnot_circuit,
# _get_reversible_tableau,
# _get_clifford_circuit,


def test_h_gadgetisation() -> None:
    circ = (
        Circuit(4)
        .T(0)
        .CX(0, 3)
        .CX(2, 1)
        .CX(3, 1)
        .T(3)
        .H(0)
        .H(1)
        .CZ(0, 3)
        .H(2)
        .CRy(0.25, 0, 3)
    )
    n_qubits_without_ancillas = circ.n_qubits
    DecomposeBoxes().apply(circ)
    ComposePhasePolyBoxes().apply(circ)
    n_internal_h_gates = get_n_internal_hadamards(circ)
    REPLACE_HADAMARDS.apply(circ)
    assert get_n_conditional_paulis(circ) == n_internal_h_gates
    assert circ.n_qubits == n_qubits_without_ancillas + n_internal_h_gates


def test_circuit_utils() -> None:
    circ = (
        Circuit(2)
        .CX(0, 1)
        .Rz(1 / 4, 1)
        .CX(0, 1)
        .Rz(-1 / 4, 1)
        .CX(0, 1)
        .Rz(0.75, 1)
        .CX(0, 1)
        .Rz(-0.75, 0)
        .Rz(-1.25, 1)
    )
    assert check_rz_angles(circ)
    circ.Rz(0.61, 0)
    assert not check_rz_angles(circ)


def test_simple_circuit() -> None:
    circ = Circuit(3)
    circ.CX(0, 1).T(1).CX(0, 1).H(0).CX(0, 2).T(2).CX(0, 2).CX(0, 1).T(1).H(0).CX(0, 1)
    draw(circ)
    ComposePhasePolyBoxes().apply(circ)
    assert circ.n_gates_of_type(OpType.H) == 2
    n_internal_h_gates = get_n_internal_hadamards(circ)
    assert n_internal_h_gates == 2
    draw(circ)
    REPLACE_HADAMARDS.apply(circ)
    assert circ.n_qubits == 5
    assert get_n_conditional_paulis(circ) == n_internal_h_gates
    draw(circ)
    assert get_n_conditional_paulis(circ) == 2
    PROPAGATE_TERMINAL_PAULI.apply(circ)
    assert get_n_conditional_paulis(circ) == 1
    draw(circ)
    # PROPAGATE_TERMINAL_PAULI.apply(circ)
    # draw(circ)

    # circ = Circuit(3)
    # circ.CX(0, 1).T(1).CX(0, 1).H(0).CX(0, 2).T(2).CX(0, 2).CX(0, 1).T(1).H(0).CX(0, 1)
    # ComposePhasePolyBoxes().apply(circ)


def test_clifford_generation() -> None:

    cnot_rz_circ = circuit_from_qasm("cnot_t_2.qasm")
    draw(cnot_rz_circ)

    qpt = QubitPauliTensor(
        qubits=cnot_rz_circ.qubits,
        paulis=[Pauli.I, Pauli.X, Pauli.I],
        coeff=1,
    )

    # Turn into a PhasePolyBox
    ComposePhasePolyBoxes().apply(cnot_rz_circ)

    phase_poly_box = cnot_rz_circ.get_commands()[0].op

    new_pauli_circ, new_s_circ = _get_circuit_fragments(
        pbox=phase_poly_box, input_pauli=qpt
    )

    # draw(new_pauli_circ)
    # draw(new_s_circ)
    # BUG - new_s_circ is not Clifford


test_clifford_generation()
