from __future__ import annotations

from pytket import Qubit
from pytket._tket.circuit import Circuit, OpType, PauliExpCommutingSetBox
from pytket.circuit import PhasePolyBox
from pytket.extensions.qiskit import qiskit_to_tk
from pytket.passes import DecomposeBoxes, ComposePhasePolyBoxes
from pytket.pauli import Pauli, QubitPauliTensor
from pytket.tableau import UnitaryTableau
from qiskit.synthesis import synth_cnot_count_full_pmh

from typing import TypeVar


from pytket.circuit.display import view_browser as draw

PAULI_DICT = {
    Pauli.I: OpType.noop,
    Pauli.X: OpType.X,
    Pauli.Y: OpType.Y,
    Pauli.Z: OpType.Z,
}


def pauli_tensor_to_circuit(pauli_tensor: QubitPauliTensor) -> Circuit:
    """Create a Circuit comprised of single qubit Pauli ops from a QubitPauliTensor."""
    pauli_circ = Circuit(len(pauli_tensor.string.to_list()))

    for qubit, pauli_op in pauli_tensor.string.map.items():
        pauli_circ.add_gate(PAULI_DICT[pauli_op], [qubit])

    return pauli_circ


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


def _get_updated_paulis(
    pauli_tensors: list[QubitPauliTensor],
    new_pauli: QubitPauliTensor,
) -> list[QubitPauliTensor]:

    new_tensors = []

    for pauli_op in pauli_tensors:
        if not pauli_op.commutes_with(new_pauli):
            new_coeff = pauli_op.coeff * (-1)
            new_tensors.append(
                QubitPauliTensor(string=pauli_op.string, coeff=new_coeff),
            )

    return new_tensors


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


def get_pauli_conjugate(
    pbox: PhasePolyBox,
    input_pauli: QubitPauliTensor,
) -> QubitPauliTensor:
    """Given a PhasePolyBox (U) and a QubitPauliTensor (P), returns P' = L P L†."""
    # Get L as a Tableau
    l_tableau = _get_reversible_tableau(pbox)

    return l_tableau.get_row_product(input_pauli)


def _get_daggered_phasepolybox(pbox: PhasePolyBox) -> PhasePolyBox:

    circuit = pbox.get_circuit()

    dg_circ = circuit.dagger()

    ComposePhasePolyBoxes().apply(dg_circ)

    return dg_circ.get_commands()[0].op


T = TypeVar("T")


def _merge_lists(list1: list[T], list2: list[T]) -> list[T]:
    assert len(list1) == len(list2)

    zipped_lists = tuple(zip(list1, list2))

    result_list = []
    for _ in range(len(zipped_lists)):
        result_list.extend(zipped_lists[_])

    return result_list


def synthesise_clifford(pbox: PhasePolyBox, input_pauli: QubitPauliTensor) -> Circuit:
    """Synthesise a Circuit implementing the end of Circuit Clifford Operator C."""
    # Get P' = L * P * L†
    new_pauli: QubitPauliTensor = get_pauli_conjugate(pbox, input_pauli)

    # Create a Circuit with the Pauli tensor P'
    result_circ: Circuit = pauli_tensor_to_circuit(new_pauli)

    # Convert U to a list of QubitPauliTensors
    s_sequence: list[QubitPauliTensor] = _parities_to_pauli_tensors(pbox)

    # Get updated S' sequence
    s_prime_sequence: list[QubitPauliTensor] = _get_updated_paulis(
        s_sequence,
        new_pauli,
    )

    # Get the PhasePolyBox implementing U†
    u_dg_box: PhasePolyBox = _get_daggered_phasepolybox(pbox)

    # U† as a list of QubitPauliTensors
    s_dg_sequence: list[QubitPauliTensor] = _parities_to_pauli_tensors(u_dg_box)

    # Get Q sequence by combining operators for S' and S† (Q = S' * S†)
    q_operator_list: list[QubitPauliTensor] = s_prime_sequence + s_dg_sequence

    # Synthesise a Phase gadget sequence for Q (Angles should be Clifford)
    q_operator_circ = _get_phase_gadget_circuit(q_operator_list)

    # Combine circuits for P' and Q
    result_circ.add_circuit(q_operator_circ, result_circ.qubits)

    return result_circ
