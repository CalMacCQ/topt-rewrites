from pytket.circuit import Circuit, OpType
from pytket.extensions.offline_display import view_browser as draw
from pytket.passes import ComposePhasePolyBoxes

from topt_rewrites.main import (
    PROPAGATE_TERMINAL_PAULI,
    REPLACE_HADAMARDS,
    check_rz_angles,
    get_n_conditional_paulis,
    get_n_internal_hadamards,
)

# def test_h_gadgetisation() -> None:
#    circ = (
#        Circuit(4)
#        .T(0)
#        .CX(0, 3)
#        .CX(2, 1)
#        .CX(3, 1)
#        .T(3)
#        .H(0)
#        .H(1)
#        .CZ(0, 3)
#        .H(2)
#        .CRy(0.25, 0, 3)
#    )
#    n_qubits_without_ancillas = circ.n_qubits
#    draw(circ)
#    DecomposeBoxes().apply(circ)
#    ComposePhasePolyBoxes().apply(circ)
#    draw(circ)
#    n_internal_h_gates = get_n_internal_hadamards(circ)
#    REPLACE_HADAMARDS.apply(circ)
#    draw(circ)
#    assert get_n_conditional_paulis(circ) == n_internal_h_gates
#    assert circ.n_qubits == n_qubits_without_ancillas + n_internal_h_gates
#


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
    # draw(circ)
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


test_simple_circuit()
