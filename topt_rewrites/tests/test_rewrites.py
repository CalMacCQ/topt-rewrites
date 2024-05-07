from pytket.circuit import Circuit, OpType, CircBox
from pytket.circuit.display import view_browser
from topt_rewrites.main import (
    replace_hadamards,
    propogate_final_pauli_x_gate,
    get_last_pauli_x_args,
)
from pytket.passes import ComposePhasePolyBoxes, DecomposeBoxes


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
    # DecomposeBoxes().apply(circ)
    # ComposePhasePolyBoxes().apply(circ)
    replace_hadamards.apply(circ)
    view_browser(circ)
    args = get_last_pauli_x_args(circ)
    print(args)
    assert circ.n_qubits == 7
    assert circ.n_bits == 3
    assert circ.n_gates_of_type(OpType.H) == 0


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
    propogate_final_pauli_x_gate(circ)


if __name__ == "__main__":
    test_hadamard_gadgetisation()
