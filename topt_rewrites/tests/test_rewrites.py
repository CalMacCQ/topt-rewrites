from pytket.circuit import Circuit, OpType, CircBox
from topt_rewrites.main import (
    replace_hadamards,
    propogate_terminal_pauli_x_gate,
    get_n_conditional_paulis,
    propogate_all_terminal_paulis,
    propogate_terminal_pauli,
)
from pytket.passes import ComposePhasePolyBoxes, DecomposeBoxes

from pytket.extensions.offline_display import view_browser as draw


def test_hadamard_gadgetisation() -> None:
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
    cb = CircBox(Circuit(2, "TEST BOX").Ry(0.29, 1).CX(1, 0))
    circ.add_gate(cb, [0, 1])
    replace_hadamards.apply(circ)
    assert circ.n_qubits == 7
    assert circ.n_bits == 3
    # assert circ.n_gates_of_type(OpType.H) == 0
    print(get_n_conditional_paulis(circ))
    draw(circ)


def test_pauli_pushing() -> None:
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
    DecomposeBoxes().apply(circ)
    ComposePhasePolyBoxes().apply(circ)
    replace_hadamards.apply(circ)
    draw(circ)
    propogate_terminal_pauli.apply(circ)
    draw(circ)
    propogate_terminal_pauli.apply(circ)
    draw(circ)


def test_simple_circuit() -> None:
    circ = Circuit(3)
    circ.CX(0, 1).T(1).CX(0, 1).H(0).CX(0, 2).T(2).CX(0, 2).CX(0, 1).T(1).H(0).CX(0, 1)
    draw(circ)
    ComposePhasePolyBoxes().apply(circ)
    assert circ.n_gates_of_type(OpType.H) == 2
    draw(circ)
    replace_hadamards.apply(circ)
    assert circ.n_qubits == 5
    draw(circ)
    assert get_n_conditional_paulis(circ) == 2
    propogate_terminal_pauli.apply(circ)
    assert get_n_conditional_paulis(circ) == 1
    draw(circ)
    propogate_terminal_pauli.apply(circ)
    draw(circ)


if __name__ == "__main__":
    test_hadamard_gadgetisation()
    test_pauli_pushing()
    test_simple_circuit()
