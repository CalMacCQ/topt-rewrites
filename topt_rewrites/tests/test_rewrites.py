from pytket import Circuit
from pytket.circuit.display import view_browser


from topt_rewrites.rewrites import replace_hadamards


def test_hadamard_gadgetisation() -> None:
    circ = Circuit(4).T(0).CX(0, 3).CX(2, 1).CX(3, 1).T(3).H(0).H(1).CZ(0, 3).H(2)
    replace_hadamards.apply(circ)
    assert circ.n_qubits == 7
    assert circ.n_bits == 3
