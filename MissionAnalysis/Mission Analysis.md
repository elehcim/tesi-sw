Appunti per la determinazione di una missione Terra-Giove usando le LCS
=======================================================================

# Assunzioni:

 * Partire dalla Terra.
 * Partire da y=0 nel sistema sinodico
   - Quindi in pratica sto assumendo un punto di partenza.
     Devo solo aggiustare la velocità (vx e vy) di partenza.
 * Scegliere l'LCS che può essere presa sfruttando al massimo l'energia della
   Terra, cioè tale per cui la differenza tra la sua energia e quella della
   terra nel sistema sinodico è minima.
   Avendo scelto y=0 la velocità della Terra è tutta lungo l'asse y (vx=0)
 * Arrivare ad un'orbita chiusa di Giove.

## Problema:

Ho delle variabili da selezionare:
 * energia
 * tempo di partenza
 * vettore velocità

## Obiettivo

Minimizzare il Delta v per raggiungere l'orbita chiusa attorno a Giove. 

### Considerazioni

 * Esiste un'energia minima istante per istante per la LCS.
   - Al di sotto di questa energia minima la LCS non è più riconoscibile.
   - Questa energia è possibile che sia collegata all'apertura del "gate"
   STUDIARE IL GATE!

[ ] - probabilmente si possono ottenere questi numeri usando le equazioni

|    t   |       e       |    Ap. L1    |    Ap. L2    |
|-------------------------------------------------------
|  0.00  |       /       |  -1.44917    |  -1.448566   |
|  0.63  |       /       |              |              |
|  0.79  |       /       |              |              |
|  1.26  | -1.44 : -1.42 |              |              |
|  1.57  | -1.48 : -1.47 |              |              |
|  2.36  | -1.54 : -1.52 |              |              |
|  3.14  | -1.60 : -1.58 |  -1.597786   |  -1.597119   |
|  3.93  | -1.59 : -1.57 |              |              |
|  4.71  | -1.52 : -1.50 |  -1.52003    |  -1.519396   |
|  5.50  | -1.43 : -1.41 |              |              |

 * Voglio partire dalla Terra perciò la LCS deve contenere la x della Terra.
   - Deve essere x_max_LCS > -0.19
   Quindi bisogna escludere i tempi e le energie per le quali questo non
   è verificato.

 * Possibilità da considerare:
   - Quanto tempo sta aperto il gate di L1 in modo da arrivare a Giove

# Inizio

 * Partenza dalla Terra *Sistema geocentrico*
 * Do un eccesso iperbolico. -> Bisogna scegliere un vettore velocità.
 * Trasformare la velocità nel sistema eliocentrico

 * Prendere l'LCS da un'altra parte uscendo dalla SOI della Terra con
 condizioni di velocità e posizione idonee


# Sistema dei Tre corpi ellittico ristretto

Propago la traiettoria

# Fine

Arrivo nella sfera di influenza di Giove

 - Energia per chiudere l'orbita nel sistema dei due corpi (energia < 0)
 - (Energia per chiudere tutti i gateway nelle condizioni peggiori di
   anomalia vera.)

# Dati

## Dimensional

Sun radius               6.96e5 km
Jupiter SOI (Sun)        48223000 km
Eart SOI (Sun)           924000 km
Jupiter-Sun distance     778547200 km
Earth orbit radius       149600000 km (=1AU)
Earth orbital speed      29.78 km/s
Jupiter orbital speed    13.07 km/s

## Non-dimensional Jupiter-Sun sinodic reference frame

mu                    9.537e-4
eccentricity          0.04839
Sun radius            8.9e-4
Jupiter SOI           0.0619
Earth SOI             0.0012
Earth orbit radius    0.1922
Earth energy          -2.80286
Earth orbital speed   2.27
Jupiter orbital speed 1.0

## Combinazioni tempo - energia

 * Apertura gate L1
        t     |    e
    -----------------
      3.14    |  -1.597786
       0      | -1.44917

 * Apertura gate L2
        t     |    e
    -----------------
      3.14    |  -1.597119
      0       |  -1.448566

 * Sparizione delle zone di Hill
        t     |    e
    -----------------
      3.14    |  -1.5770
      0       |  -1.4303
