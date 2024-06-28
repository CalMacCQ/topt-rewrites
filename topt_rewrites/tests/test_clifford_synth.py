from __future__ import annotations

from glob import glob

from pytket.circuit import PhasePolyBox
from pytket.circuit.display import view_browser as draw
from pytket.passes import ComposePhasePolyBoxes
from pytket.predicates import CliffordCircuitPredicate
from pytket.qasm import circuit_from_qasm

from topt_rewrites.clifford import synthesise_clifford
from topt_rewrites.utils import tensor_from_x_index


def test_clifford_synthesis(
    x_index: int,
    draw_circuit: bool,
    circuit_list: list[str],
) -> None:
    """Tests Clifford synthesis method."""
    for circuit_file in circuit_list:
        cnot_rz_circ = circuit_from_qasm(circuit_file)

        print("=========================================")
        print(f"Now Testing Circuit: {circuit_file}")

        if draw_circuit:
            draw(cnot_rz_circ)

        # Turn into a PhasePolyBox
        ComposePhasePolyBoxes().apply(cnot_rz_circ)

        phase_poly_box: PhasePolyBox = cnot_rz_circ.get_commands()[0].op

        qpt = tensor_from_x_index(x_index=x_index, n_qubits=cnot_rz_circ.n_qubits)

        result_circ = synthesise_clifford(pbox=phase_poly_box, input_pauli=qpt)
        print(
            f"Is synthesised circuit for {circuit_file} Clifford?",
            CliffordCircuitPredicate().verify(result_circ),
        )

        if draw_circuit:
            draw(result_circ)


if __name__ == "__main__":
    CIRCUIT_LIST = glob("qasm/*.qasm")
    X_INDEX = 1
    DRAW_CIRCUIT = False
    test_clifford_synthesis(
        x_index=X_INDEX,
        draw_circuit=DRAW_CIRCUIT,
        circuit_list=CIRCUIT_LIST,
    )


# circ1 = Circuit(2).CX(1, 0).T(0).CX(1, 0).X(1).CX(1, 0).Tdg(0).CX(1, 0)
#
# circ2 = Circuit(2).X(1).CX(1, 0).Sdg(0).CX(1, 0)
#
# u1 = circ1.get_unitary()
#
# u2 = circ2.get_unitary()
#
# from pytket.utils import compare_unitaries
#
# print(compare_unitaries(u1, u2))
