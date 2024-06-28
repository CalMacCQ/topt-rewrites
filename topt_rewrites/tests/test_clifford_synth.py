from __future__ import annotations

from glob import glob

from pytket.circuit import PhasePolyBox, Circuit
from pytket.circuit.display import view_browser as draw
from pytket.passes import ComposePhasePolyBoxes
from pytket.pauli import QubitPauliTensor
from pytket.predicates import CliffordCircuitPredicate
from pytket.qasm import circuit_from_qasm
from pytket.utils import compare_unitaries

from topt_rewrites.clifford import synthesise_clifford, pauli_tensor_to_circuit
from topt_rewrites.utils import tensor_from_x_index


def get_test_circuit(p_box: PhasePolyBox, pauli: QubitPauliTensor) -> Circuit:
    circ = Circuit(p_box.n_qubits)
    u_circ = p_box.get_circuit()
    circ.append(u_circ)
    pauli_circ = pauli_tensor_to_circuit(pauli)
    circ.append(pauli_circ)
    u_dg_circ = u_circ.dagger()
    circ.append(u_dg_circ)
    return circ


def test_clifford_synthesis(
    x_index: int,
    circuit_list: list[str],
    draw_circuit: bool = False,
    check_unitary: bool = False,
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

        if check_unitary:
            test_circuit = get_test_circuit(phase_poly_box, qpt)
            draw(test_circuit)
            test_unitary = test_circuit.get_unitary()
            result_unitary = result_circ.get_unitary()
            print(
                f"Do the unitaries match? {compare_unitaries(test_unitary, result_unitary)}",
            )


if __name__ == "__main__":
    # CIRCUIT_LIST = glob("qasm/*.qasm")
    CIRCUIT_LIST = ["qasm/cnot_t_0.qasm"]
    X_INDEX = 1
    DRAW_CIRCUITS = False
    CHECK_UNITARY = True
    test_clifford_synthesis(
        x_index=X_INDEX,
        draw_circuit=DRAW_CIRCUITS,
        circuit_list=CIRCUIT_LIST,
        check_unitary=CHECK_UNITARY,
    )
