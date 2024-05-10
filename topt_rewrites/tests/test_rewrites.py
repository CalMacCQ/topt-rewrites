from pytket.circuit import Circuit, OpType, CircBox
from topt_rewrites.main import (
    REPLACE_HADAMARDS,
    get_n_conditional_paulis,
    # PROPOGATE_TERMINAL_PAULI,
    get_n_internal_hadamards,
)
from pytket.passes import ComposePhasePolyBoxes, DecomposeBoxes

# from pytket.extensions.offline_display import view_browser as draw


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
    # draw(circ)
    DecomposeBoxes().apply(circ)
    ComposePhasePolyBoxes().apply(circ)
    # draw(circ)
    n_internal_h_gates = get_n_internal_hadamards(circ)
    REPLACE_HADAMARDS.apply(circ)
    assert get_n_conditional_paulis(circ) == n_internal_h_gates
    assert circ.n_qubits == n_qubits_without_ancillas + n_internal_h_gates


def test_simple_circuit() -> None:
    circ = Circuit(3)
    circ.CX(0, 1).T(1).CX(0, 1).H(0).CX(0, 2).T(2).CX(0, 2).CX(0, 1).T(1).H(0).CX(0, 1)
    # draw(circ)
    ComposePhasePolyBoxes().apply(circ)
    assert circ.n_gates_of_type(OpType.H) == 2
    n_internal_h_gates = get_n_internal_hadamards(circ)
    assert n_internal_h_gates == 2
    # draw(circ)
    REPLACE_HADAMARDS.apply(circ)
    assert circ.n_qubits == 5
    assert get_n_conditional_paulis(circ) == n_internal_h_gates
    # draw(circ)
    # assert get_n_conditional_paulis(circ) == 2
    # PROPOGATE_TERMINAL_PAULI.apply(circ)
    # assert get_n_conditional_paulis(circ) == 1
    # draw(circ)
    # PROPOGATE_TERMINAL_PAULI.apply(circ)
    # draw(circ)


if __name__ == "__main__":
    test_hadamard_gadgetisation()
    test_pauli_pushing()
    test_simple_circuit()
