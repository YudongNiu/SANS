# SCANS: Efficient Densest Subgraph Discovery over Relational Graphs without Materialization

Compile
-------
Before compile the program, set the key parameters for our program in file "param.h".

**There are 5 parameters:**
* K: the size of the KMV sketch
* Kmin: reconstruction threshold
* L: the number of KMV sketch constructed initially

Then, compile the code for SCANS by executing the following command on linux:

```sh
g++ -O3 main.py -o scans
```

Running code
-------

To run the code for densest subgraph discovery queries, execute the following command on linux:

```sh
./scans dataset method
```

**There are 5 parameters:**
* dataset: name of the KG
* method: the algorithm to run

For example, the following command execute the ScansE method in each relational graph induced by candidate meta-paths in imdb.

```sh
./scans imdb scanse
```

Input Files
-----------
**The program scans requires 4 input files:**
* node.dat stores nodes in KG.
* link.dat stores edges in KG.
* meta.dat stores the number of nodes in KG.
* dataset-cod-global-rules.dat stores the meta-paths mined by AnyBURL
