# pgraph-tascel

Sequence homology detection is central to a number of bioinformatics applications including genome sequencing and protein family characterization. Given millions of sequences, the goal is to identify all pairs of sequences that are highly similar (or “homologous”) on the basis of alignment criteria. While there are optimal alignment algorithms to compute pairwise homology, their deployment for large-scale is currently not feasible; instead, heuristic methods are used at the expense of quality. Here, we present a parallel implementation for conducting optimal homology detection on distributed memory supercomputers. Our approach uses a combination of techniques from asynchronous load balancing (viz. work stealing, dynamic task counters), data replication, and exact-matching filters to achieve homology detection at scale.

## What This Softare Does

This software performs pairwise sequence alignment while balancing the load across MPI ranks.  The load imbalance arises from the cost model of pairwise sequence alignment depending on the length of the two input sequences.  The load is further imbalanced because a preprocessing step selects for alignment only the sequences containing an exact match of a given length.  The original filter used a suffix tree data structure.  However, this was less efficient than using a suffix array data structure, so the latest code was eventually converted to use the latter.  There were many load balancing strategies attempted, notably work stealing and distributed task counters (see reference below).  The work stealing approach, though promising, was never released publicly.  The remaining code using a distributed task counter is described in the next section.

Note that this software relies on the [parasail](https://github.com/jeffdaily/parasail) code for performing sequence alignments.  In the end, this means pgraph-tascel really only adds a suffix array filter and task counter to the pipeline of sequence alignment.

## How to Use the MPI-only Code

The code that is free from unusable dependencies is the `align_parted_nxtval` code in the apps directory.  It relies on a distributed task counter that lives on MPI rank 0 in a separate thread of execution.  This requires an MPI implementation that provides MPI_THREAD_MULTIPLE support.  The input parameters are specified using a YAML file.  It is highly recommended to run a single MPI rank per node, otherwise known as MPI+X where X in this case will be OpenMP.

The input FASTA file is broadcast to all MPI ranks.  The all-to-all sequence alignment is broken up into tiles, and each tile represents a task in the task counter.  There are two types of tiles, those representing sequence sets that are compared with themselves and those representing a sequence set that is compared against a different sequence set.  For each tile, a suffix array is constructed for the sequences represented by the tile.  The sequence pairs that are not filtered out by the suffix array are then aligned using an OpenMP loop and the parasail software.

The output of the softare is a series of graph edges where highly similar sequences are the vertexes of the edge.  The graph can then be postprocessed and run through the grappolo code https://github.com/luhowardmark/GrappoloTK.

The parasail software referred to above contains an application, the parasail_aligner, that performs the same steps as the align_parted_nxtval code but does not use MPI (it runs on a single workstation) and treats the input as a single tile.  In other words, the parasail_aligner constructs the suffix array filter on the entire set of sequences, performs the alignments, and outputs a graph directly.

## How to Use the Old Tascel-based Code

Much of this software depends on a software library that was never released publicly, 'tascel 2.0'.  The predecessor to tascel 2.0 was released as version 0.1 with the Global Arrays project - https://github.com/GlobalArrays/tascel - which in turn depends on the ARMCI library.  The ARMCI library builds and installs as part of the Global Arrays library - https://github.com/GlobalArrays/ga.  There is a chance that a motivated user might be able to adapt the pgraph-tascel code back to using only the APIs from tascel 0.1.

## Citing pgraph-tascel

If needed, please cite the following paper.

Daily JA, A Kalyanaraman, S Krishnamoorthy, and A Vishnu. 2015.
"A Work Stealing Based Approach for Enabling Scalable Optimal Sequence Homology Detection."
Journal of Parallel and Distributed Computing 79-80:132-142.
doi:10.1016/j.jpdc.2014.08.009

http://dx.doi.org/10.1016/j.jpdc.2014.08.009

## License: Battelle BSD-style

Copyright (c) 2015, Battelle Memorial Institute

1.  Battelle Memorial Institute (hereinafter Battelle) hereby grants
    permission to any person or entity lawfully obtaining a copy of this
    software and associated documentation files (hereinafter “the
    Software”) to redistribute and use the Software in source and binary
    forms, with or without modification.  Such person or entity may use,
    copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and may permit others to do so, subject to
    the following conditions:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimers.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.

    - Other than as used herein, neither the name Battelle Memorial
      Institute or Battelle may be used in any form whatsoever without
      the express written consent of Battelle.

    - Redistributions of the software in any form, and publications
      based on work performed using the software should include the
      following citation as a reference:

      Daily JA, A Kalyanaraman, S Krishnamoorthy, and A Vishnu. 2015.
      "A Work Stealing Based Approach for Enabling Scalable Optimal Sequence Homology Detection."
      Journal of Parallel and Distributed Computing 79-80:132-142.
      doi:10.1016/j.jpdc.2014.08.009

2.  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL BATTELLE
    OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
    USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
    OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
    SUCH DAMAGE.
