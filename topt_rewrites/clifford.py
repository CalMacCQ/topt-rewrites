from __future__ import annotations

from pytket._tket.circuit import Circuit, OpType, PauliExpCommutingSetBox

from pytket import Qubit
from pytket.circuit import PhasePolyBox
from pytket.extensions.qiskit import qiskit_to_tk
from pytket.passes import DecomposeBoxes
from pytket.pauli import Pauli, QubitPauliTensor
from pytket.tableau import UnitaryTableau
from qiskit.synthesis import synth_cnot_count_full_pmh

PAULI_DICT = {
    Pauli.I: OpType.noop,
    Pauli.X: OpType.X,
    Pauli.Y: OpType.Y,
    Pauli.Z: OpType.Z,
}


def _get_reversible_tableau(pbox: PhasePolyBox) -> UnitaryTableau:
    # cheat by synthesising the CNOT circuit with qiskit and converting
    qc = synth_cnot_count_full_pmh(pbox.linear_transformation, section_size=2)
    qc2 = qc.reverse_bits()  # correct for endianness
    tkc_cnot = qiskit_to_tk(qc2)
    return UnitaryTableau(tkc_cnot)


def _parities_to_pauli_tensors(pbox: PhasePolyBox) -> list[QubitPauliTensor]:
    phase_poly = pbox.phase_polynomial
    tensor_list = []
    for parity, phase in phase_poly.items():
        pauli_list = []
        qubit_list = []
        for count, boolean in enumerate(parity):
            qubit_list.append(Qubit(count))
            if boolean:
                pauli_list.append(Pauli.Z)
            else:
                pauli_list.append(Pauli.I)

        pauli_tensor = QubitPauliTensor(
            qubits=qubit_list,
            paulis=pauli_list,
            coeff=phase,
        )
        tensor_list.append(pauli_tensor)

    return tensor_list


def get_circuit_fragments(
    pbox: PhasePolyBox,
    input_pauli: QubitPauliTensor,
) -> tuple[Circuit, Circuit]:
    # p_box = S L
    # S: Sequence of {Z, I} Phase gadgets
    # L: Linear reversible transformation

    # Compute Tableau form of L.
    l_tableau: UnitaryTableau = _get_reversible_tableau(pbox)

    # Updated Pauli Tensor P'. P' = L P L†.
    new_pauli: QubitPauliTensor = l_tableau.get_row_product(input_pauli)

    # Create new circuit implementing P'
    pauli_prime_circ = Circuit(pbox.n_qubits)
    for qubit, pauli_op in new_pauli.string.map.items():
        pauli_prime_circ.add_gate(PAULI_DICT[pauli_op], [qubit])

    # Convert (True, False) terms to Pauli.Z, Pauli.I in a QubitPauliTensor.
    pauli_tensors = _parities_to_pauli_tensors(pbox)

    # Get updated phase gadget Sequence S' and synthesise circuit.
    s_prime = _get_updated_paulis(pauli_tensors, new_pauli)
    s_prime_circ = _get_phase_gadget_circuit(s_prime)

    return pauli_prime_circ, s_prime_circ


# def get_clifford_circuit()


def _get_updated_paulis(
    pauli_tensors: list[QubitPauliTensor],
    new_pauli: QubitPauliTensor,
) -> list[QubitPauliTensor]:

    for pauli_op in pauli_tensors:
        if not pauli_op.commutes_with(new_pauli):
            pauli_op.coeff *= 2

    return pauli_tensors


def _get_phase_gadget_circuit(pauli_tensors: list[QubitPauliTensor]) -> Circuit:

    pauli_ops = []
    for tensor in pauli_tensors:
        pauli_list = list(tensor.string.map.values())
        pair = (pauli_list, tensor.coeff.real)
        pauli_ops.append(pair)

    pauli_gadgets_sequence = PauliExpCommutingSetBox(pauli_ops)
    n_qubits = pauli_gadgets_sequence.n_qubits

    pauli_gadget_circ = Circuit(n_qubits).add_gate(
        pauli_gadgets_sequence,
        list(range(n_qubits)),
    )

    DecomposeBoxes().apply(pauli_gadget_circ)

    return pauli_gadget_circ
