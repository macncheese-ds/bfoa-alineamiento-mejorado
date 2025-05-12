import copy
import math
import random
from multiprocessing import Manager, Pool
from evaluadorBlosum import evaluadorBlosum
import numpy
from fastaReader import fastaReader
import concurrent.futures


class bacteria:
    def __init__(self, numBacterias, initial_step=1):
        manager = Manager()
        self.blosumScore = manager.list([0] * numBacterias)
        self.tablaAtract = manager.list([0] * numBacterias)
        self.tablaRepel = manager.list([0] * numBacterias)
        self.tablaInteraction = manager.list([0] * numBacterias)
        self.tablaFitness = manager.list([0] * numBacterias)
        self.granListaPares = manager.list([[] for _ in range(numBacterias)])
        self.NFE = manager.list([0] * numBacterias)
        self.step_sizes = manager.list([initial_step] * numBacterias)

    def resetListas(self, numBacterias, initial_step=1):
        manager = Manager()
        self.blosumScore = manager.list([0] * numBacterias)
        self.tablaAtract = manager.list([0] * numBacterias)
        self.tablaRepel = manager.list([0] * numBacterias)
        self.tablaInteraction = manager.list([0] * numBacterias)
        self.tablaFitness = manager.list([0] * numBacterias)
        self.granListaPares = manager.list([[] for _ in range(numBacterias)])
        self.NFE = manager.list([0] * numBacterias)
        self.step_sizes = manager.list([initial_step] * numBacterias)

    def _mueve_gaps(self, bacterSeq, numGaps):
        """
        Inserta numGaps gaps ('-') en posiciones aleatorias de las secuencias de bacterSeq.
        bacterSeq es una lista de listas (cada lista es una secuencia).
        """
        seq_mod = copy.deepcopy(bacterSeq)
        for _ in range(numGaps):
            # elige una de las secuencias al azar
            seq_idx = random.randrange(len(seq_mod))
            # elige posición de inserción
            pos = random.randrange(len(seq_mod[seq_idx]) + 1)
            seq_mod[seq_idx] = seq_mod[seq_idx][:pos] + ['-'] + seq_mod[seq_idx][pos:]
        return seq_mod

    def tumbo(self, numSec, poblacion, numGaps):
        for i in range(len(poblacion)):
            bacterTmp = list(poblacion[i])
            # delegar inserción de gaps
            bacterTmp = self._mueve_gaps(bacterTmp, numGaps)
            poblacion[i] = tuple(bacterTmp)

    def nado(self, numSec, poblacion, step_size):
        for i in range(len(poblacion)):
            pos_actual = list(poblacion[i])
            fitness_actual = self.evaluaFila(pos_actual, i)
            for swim in range(3):
                # usar paso adaptativo si existe
                gaps = step_size
                if hasattr(self, 'step_sizes'):
                    gaps = self.step_sizes[i]
                nuevo = self._mueve_gaps(pos_actual, int(gaps))
                fit_nuevo = self.evaluaFila(nuevo, i)
                if fit_nuevo > fitness_actual:
                    poblacion[i] = tuple(nuevo)
                    pos_actual = nuevo
                    fitness_actual = fit_nuevo
                else:
                    break

    def ajusta_paso(self, iter, max_iter, decay=0.1):
        """
        Reduce el tamaño de paso de cada bacteria linealmente según iteración.
        """
        for i in range(len(self.step_sizes)):
            self.step_sizes[i] *= (1 - decay * (iter / max_iter))

    def elitismo(self, poblacion, tasa=0.1):
        """
        Conserva el top tasa% de bacterias y reemplaza a las peores.
        """
        # calcula índices ordenados por fitness
        indices = sorted(range(len(self.tablaFitness)), key=lambda i: self.tablaFitness[i], reverse=True)
        num_elite = max(1, int(len(poblacion) * tasa))
        elites = [poblacion[i] for i in indices[:num_elite]]
        # reemplaza peores
        for j, idx in enumerate(indices[-num_elite:]):
            poblacion[idx] = copy.deepcopy(elites[j])

    def mutacion(self, poblacion, prob=0.05, numGaps=1):
        """
        Aplica mutación aleatoria con probabilidad prob.
        """
        for i in range(len(poblacion)):
            if random.random() < prob:
                seq = list(poblacion[i])
                seq = self._mueve_gaps(seq, numGaps)
                poblacion[i] = tuple(seq)

    # Métodos originales (sin cambios) abajo...

    def cuadra(self, numSec, poblacion):
        for i in range(len(poblacion)):
            bacterTmp = list(poblacion[i])
            bacterTmp = bacterTmp[:numSec]
            maxLen = 0
            for j in range(numSec):
                if len(bacterTmp[j]) > maxLen:
                    maxLen = len(bacterTmp[j])
                    for t in range(numSec):
                        gap_count = maxLen - len(bacterTmp[t])
                        if gap_count > 0:
                            bacterTmp[t].extend(["-"] * gap_count)
                            poblacion[i] = tuple(bacterTmp)

    def creaGranListaPares(self, poblacion):
        for i in range(len(poblacion)):
            pares = []
            bacterTmp = list(poblacion[i])
            for j in range(len(bacterTmp)):
                column = self.getColumn(bacterTmp, j)
                pares += self.obtener_pares_unicos(column)
            self.granListaPares[i] = pares

    def evaluaFila(self, fila, num):
        evaluador = evaluadorBlosum()
        score = 0
        for par in fila:
            score += evaluador.getScore(par[0], par[1])
        self.blosumScore[num] = score
        return score

    def evaluaBlosum(self):
        with Pool() as pool:
            args = [(copy.deepcopy(self.granListaPares[i]), i) for i in range(len(self.granListaPares))]
            pool.starmap(self.evaluaFila, args)

    def getColumn(self, bacterTmp, colNum):
        return [seq[colNum] for seq in bacterTmp]

    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                pares_unicos.add(tuple(sorted([columna[i], columna[j]])))
        return list(pares_unicos)

    def compute_diff(self, args):
        indexBacteria, otherBlosumScore, selfScores, d, w = args
        diff = (selfScores[indexBacteria] - otherBlosumScore) ** 2.0
        self.NFE[indexBacteria] += 1
        return d * numpy.exp(w * diff)

    def compute_cell_interaction(self, indexBacteria, d, w, atracTrue):
        with Pool() as pool:
            args = [(indexBacteria, other, self.blosumScore, d, w) for other in self.blosumScore]
            results = pool.map(self.compute_diff, args)
            pool.close()
            pool.join()
        total = sum(results)
        if atracTrue:
            self.tablaAtract[indexBacteria] = total
        else:
            self.tablaRepel[indexBacteria] = total
    
    def creaTablaAtract(self, poblacion, d, w):
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria, d, w, atracTrue=True)

    def creaTablaRepel(self, poblacion, d, w):
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria, d, w, atracTrue=False)


    def creaTablaAtractRepel(self, poblacion, dAttr, wAttr, dRepel, wRepel):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.submit(self.creaTablaAtract, poblacion, dAttr, wAttr)
            executor.submit(self.creaTablaRepel, poblacion, dRepel, wRepel)

    def creaTablaInteraction(self):
        for i in range(len(self.tablaAtract)):
            self.tablaInteraction[i] = self.tablaAtract[i] + self.tablaRepel[i]

    def creaTablaFitness(self):
        for i in range(len(self.tablaInteraction)):
            self.tablaFitness[i] = self.blosumScore[i] + self.tablaInteraction[i]

    def getNFE(self):
        return sum(self.NFE)

    def obtieneBest(self, globalNFE):
        bestIdx = max(range(len(self.tablaFitness)), key=lambda i: self.tablaFitness[i])
        print(f"Best: {bestIdx} Fitness: {self.tablaFitness[bestIdx]}")
        return bestIdx, self.tablaFitness[bestIdx]

    def replaceWorst(self, poblacion, best):
        worst = min(range(len(self.tablaFitness)), key=lambda i: self.tablaFitness[i])
        poblacion[worst] = copy.deepcopy(poblacion[best])
