# Una semplice iterazione di Jacobi

***

## Programmazione Concorrente, Parallela e su Cloud
### Università degli Studi di Salerno
#### *Anno Accademico 2017/2018*

**Docenti:** _Vittorio Scarano_ , _Carmine Spagnuolo_  
**Studente:** _Gianmarco Russo_  

---

## Problem Statement
Si vuole implementare una semplice iterazione di Jacobi per l'approssimazione di soluzioni a sistemi lineari di equazioni. Le equazioni saranno rappresentate tramite una matrice contenente i coefficienti delle incognite di ogni equazione e i termini noti.  
Il seguente codice
```
while (not converged) {
  for (i,j)
    xnew[i][j] = (x[i+1][j] + x[i-1][j] + x[i][j+1] + x[i][j-1])/4;
  for (i,j)
    x[i][j] = xnew[i][j];
  }
  ```
calcola un'approssimazione per la soluzione di un sistema di equazioni. Per questo caso di studio i valori di "bordo" sono fissi e il calcolo viene effettuato solo sui valori interni. In pratica, ciò significa che in una matrice nxm, i valori
```
x[0][j]
x[n-1][j]
x[i][0]
x[i][m-1]
```
con  0&lt;i&lt;n e 0&lt;j&lt;m, restano immutati.  
Si vuole computare  questa approssimazione in parallelo.  
Per il test di convergenza si calcola:
```
diffnorm = 0;
for (i,j)
    diffnorm += (xnew[i][j] - x[i][j]) * (xnew[i][j] - x[i][j]);
diffnorm = sqrt(diffnorm);
```
quando diffnorm è minore di 1.0e-2, si considera l'esecuzione convergente, altrimenti la si blocca se si raggiungono 100 iterazioni.  

## Soluzione e implementazione 
La soluzione proposta considera matrici qualsiasi **n x m** e un qualsiasi numero di processori **p** con **n,m,p ∈ ℤ+**. Sono state usate funzioni collettive, come richiesto, quali _Reduce_ e _Broadcast_, con l'aggiunta di funzioni _Send_ e _Receive_ bloccanti, per comunicazioni tra coppie di processori "contigui".

Per parallelizzare il problema, la dimensione della matrice è conosciuta da tutti i processori, i quali, con la seguente formula, 
```
r = n%p;
size = my_rank<r?n/p+1:n/p;
```
calcolano la quantità di righe della matrice di cui si devono occupare.  
Successivamente ognuno alloca le proprie sottomatrice in questo modo:
```
x = malloc((size+2)*m*sizeof(float));
xnew = malloc((size+2)*m*sizeof(float));
```
ognuna con 2 righe aggiuntive dove inserire le righe ricevute dal processore precedente e da quello successivo (in ordine di rango). Le sottomatrici sono in realtà dichiarate come degli array per poter essere inviate in un'unica send nella fase finale. 
Esse sono poi riempite in modo pseudocasuale. Si utilizza un offset
```
offset = my_rank<r?my_rank*size:my_rank*size+r;
```
in modo che ogni processore esegua lo stesso ciclo per generare numeri tramite la rand, ma memorizzi nella sua sottomatrice la sequenza che effettivamente avrebbe occupato quelle righe nel caso di una matrice unica su programma sequenziale (righe 56-67).
Per collezionare i risultati finali il processore con rango 0 alloca una matrice **n x m** nella quale inserire le sottomatrici che tutti gli altri processori inviano alla fine della computazione. 

Il ciclo nel quale viene eseguito il vero e proprio calcolo si può idealmente suddividere in 3 fasi. Nella prima i processori inviano la prima riga al precedente (tranne il primo) e l'ultima riga al sussessivo (tranne l'ultimo). Ricevute le righe dei "vicini", nella seconda fase, ognuno inizia la computazione sulla propria sottomatrice, calcola il diffnorm locale e copia i nuovi valori da **xnew** a **x**. Nell'ultima fase il master con una reduce calcola il diffnorm totale e ne fa la broadcast a tutti gli altri che usano questo valore per controllare la condizione del ciclo.  
Alla fine del ciclo, il master riceve le sottomatrici dagli altri processori e, insieme alla sua, le unisce nella matrice **xFinal**, che viene poi stampata a video seguita dal numero di passi effettuati, dal valore di diffnorm e dal tempo di esecuzione impiegato. 
Per la ricezione delle sottomatrici è stata utilizzata comunicazione di tipo point-to-point. Inizialmente ogni processore inviava con una prima _Send_ la dimensione della sua sottomatrice e con una seconda la sottomatrice stessa. Successivamente la prima _Send_ è stata eliminata perchè il master può calcolare da solo la dimensione della sottomatrice del processore che invia conoscendo rango, dimensione e resto. Si è preferita la comunicazione point-to-point anche per fare in modo che il master riceva le sottomatrici in ordine, così che la matrice **xFinal** rispecchi esattamente la matrice completa finale.  

### Testing
I test sono stati effettuati su istanze m4.xlarge (macchine a 4 core) di Amazon Web Service. Sono state utilizzate 8 Istanze EC2 m4.xlarge, quindi un massimo di 32 processori. Tutti i tempi raccolti tengono conto del tempo di esecuzione dall'inizio del ciclo di calcolo fino alla comunicazione dei risultati parziali al termine del ciclo. I tempi sono stati misurati in millisecondi usando la funzione _gettimeofday()_ in questo modo
```
struct timeval stop, start;

gettimeofday(&start, NULL);
//codice di calcolo
gettimeofday(&stop, NULL);
double msec = ((stop.tv_sec - start.tv_sec) * 1000.0)
    + ((stop.tv_usec - start.tv_usec) / 1000.0);
printf("time: %f ms\n",msec);
```
Come grandezza di input è stato scelto il numero di righe della matrice visto che la divisione tra i processori viene fatta, appunto, per righe. Il numero di colonne nei test è sempre lo stesso **m=1000**.

## Strong Scalabillity
Per la srong scalability è stata usata una matrice con **n=50000**.
L'esecuzione del probramma sequenziale ha impiegato **88671.449 ms**. Di seguito le esecuzioni con MPI: 

| N. processori | Tempo (millisecondi) |
|--------|--------|
| 2 |  42306,472 |
| 4 |  21838,090 |
| 8 |  11372,334 |
| 12 | 8096,803 |
| 16 | 6235,488 |
| 20 | 8940,354 |
| 24 | 7601,191 |
| 28 | 6330,158 |
| 32 | 5996,831 |

Di seguito il grafico corrispondente: 

![StrongScalability](/StrongScalability.png "Strong Scalability")

Si può notare un leggero aumento di tempo dall'utilizzo di 20 processori in poi dovuto probabilmente all'overhea di comunicazione. 

## Weak Scalability
Per la weak scalability la dimensione dell'input deve crescere proporzionalmente al numero di processori. Si è scelto quindi di definire il numero di righe della matrice in funzione di p, in questo modo:  **n=2000\*p** dove **p** è il numero di processori impiegati.
L'esecuzione del probramma sequenziale ha impiegato **3577.011 ms** (2000 righe). Di seguito le esecuzioni con MPI: 
| N. processori | Tempo (millisecondi) |
|--------|--------|
| 2  |  3401,783 |
| 4  |  3533,774 |
| 8  |  3758,116 |
| 12  | 3892,840 |
| 16  | 4058,191 |
| 20  | 7028,294 |
| 24  | 7450,617 |
| 28  | 7126,919 |
| 32  | 7224,257 |

Di seguito il grafico corrispondente: 

![WeakScalability](/WeakScalability.png "Weak Scalability")

Si può notare un leggero aumento di tempo dall'utilizzo di 20 processori, come nei test di strong scalability, dovuto probabilmente all'overhea di comunicazione. 

## Fattori di scalabilità

Per i test di Strong Scalability i fattori di scalabilità sono stati calcolati tramite la formula
```
t1/(p*tp)*100%
```
Per i test di Weak Scalability i fattori di scalabilità sono stati calcolati tramite la formula
```
(t1/tp)*100%
```
Con:
- t1: tempo di esecuzione utilizzando 1 processore
- p: numero di processori utilizzati 
- tp: tempo di esecuzione utilizzando N processori

I fattori di scalabilità sono riportati nella tabella (approssimati a tre cifre dopo la virgola):

||2|4|8|12|16|20|24|28|32|
|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
|Strong Scalability| 1,048 | 1,015 | 0,975 | 0,913 | 0,889 | 0,496 | 0,486 | 0,500 | 0,462 |
|Weak Scalability| 1,051 | 1,012 | 0,952 | 0,919 | 0,881 | 0,509 | 0,480 | 0,502 | 0,494 |

---

##### Project Template
* PCPCGianmarcoRusso.c: file C con codice per esecuzione con MPI
* PCPCSeq.c: file C con codice per esecuzione sequenziale
