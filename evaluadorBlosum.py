import blosum as bl
import numpy as np

class evaluadorBlosum:
    def __init__(self):
        # Carga la matriz BLOSUM62
        matrix = bl.BLOSUM(62)
        # Obtener alfabeto de la matriz (lista de aminoácidos)
        try:
            aa_list = matrix.alphabet
        except AttributeError:
            # Si la biblioteca usa keys() para el diccionario
            aa_list = list(matrix.keys())
        self.aa_list = aa_list
        # Mapear aminoácidos a índices
        self.aa_to_idx = {aa: idx for idx, aa in enumerate(self.aa_list)}
        # Construir matriz numpy para scores
        size = len(self.aa_list)
        self.blosum_mat = np.zeros((size, size), dtype=int)
        for i, aa1 in enumerate(self.aa_list):
            for j, aa2 in enumerate(self.aa_list):
                self.blosum_mat[i, j] = matrix[aa1][aa2]
        # Penalización por gap
        self.gap_penalty = -8

    def showMatrix(self):
        """Imprime la matriz BLOSUM vectorizada"""
        print(self.blosum_mat)

    def getScore(self, a, b=None):
        """
        Si se pasa un par de caracteres (a, b) devuelve el score simple.
        Si se pasa una lista de pares [(A1,B1), (A2,B2), ...], devuelve la suma vectorizada de scores.
        """
        # Vectorizado: lista de pares
        if isinstance(a, (list, tuple, np.ndarray)) and b is None:
            pairs = list(a)
            # Construir arrays de índices
            rows = np.array([self.aa_to_idx.get(x, -1) for x, _ in pairs])
            cols = np.array([self.aa_to_idx.get(y, -1) for _, y in pairs])
            # Aplicar penalización de gaps cuando índice == -1
            mask_gap = (rows < 0) | (cols < 0)
            scores = np.where(mask_gap, self.gap_penalty, self.blosum_mat[rows, cols])
            return int(scores.sum())
        # Caso individual
        A, B = a, b
        if A == '-' or B == '-':
            return self.gap_penalty
        idxA = self.aa_to_idx.get(A)
        idxB = self.aa_to_idx.get(B)
        if idxA is None or idxB is None:
            # Aminoácido desconocido
            return 0
        return int(self.blosum_mat[idxA, idxB])

    # Método opcional para evaluar múltiples filas de pares
    def evalua_pares(self, lista_filas):
        """
        Recibe una lista de listas de pares y devuelve una lista de scores.
        { [(A1,B1),(A2,B2)], [(C1,D1),(C2,D2)], ... }
        """
        resultados = []
        for pares in lista_filas:
            resultados.append(self.getScore(pares))
        return resultados
