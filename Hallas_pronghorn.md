# These are unix commands commonly used in Genomic analyses and HPC systems

## Screen
Screen or GNU Screen is a terminal multiplexer. In other words, it means that you can start a screen session and then open any number of windows (virtual terminals) inside that session. Processes running in Screen will continue to run when their window is not visible even if you get disconnected.
___

* To create a named session
  
```
screen -S session_name
```

* Detach
  
```
Ctrl + a + d
```
___
___

## SLURM + HPC
Slurm is an open source, fault-tolerant, and highly scalable cluster management and job scheduling system for large and small Linux clusters. Slurm requires no kernel modifications for its operation and is relatively self-contained. As a cluster workload manager, Slurm has three key functions. First, it allocates exclusive and/or non-exclusive access to resources (compute nodes) to users for some duration of time so they can perform work. Second, it provides a framework for starting, executing, and monitoring work (normally a parallel job) on the set of allocated nodes. Finally, it arbitrates contention for resources by managing a queue of pending work.  Web resources: [Slurm website](https://slurm.schedmd.com/overview.html) and [quick cheat sheet](https://slurm.schedmd.com/pdfs/summary.pdf)

*IMPORTANT: always know how much memory and CPUs are avaiable for HPC in order to optimize analyses.* 

**The BioNRES node 68 has:**
* 192 GB physical memory. 
  * However, 24 GB are reserved for processes, which means there is roughly 5GB per core.
  * Allocation of memory is important for BioNRES becuase if a user is running a job with 16 cores but uses all the memory than the other 16 cores cannot be used by another user.
* 32 cores / 64 threads
  * We can use threading on BioNRES but can't use MPI because the account only allows use of one node.  This is dependant upon how you allocate resources. (i.e. `--ntasks` and `--cpus-per-task`).  In addition, `--hint=compute_bound` allows only one thread per cpu.  With this option enabled only 32 threads can be activated.  To run more than 32 threads don't include `--hint=compute_bound`.  The benefit of `--hint=compute_bound` is that it doesn't allow another job to use a thread on a CPU currently being used, which can slow an analysis.  Or at least that is what I understand from Sebastian.  

  **It is always best to check allocation of resources with Sebastian before submitting a job.**

    *   Thread Example:
        ```
        #SBATCH --ntasks 1
        #SBATCH --cpus-per-task 32
        ```
    *   MPI Example:
         ```
        #SBATCH --ntasks 64
        #SBATCH --cpus-per-task 1
        ```

*Example of Slurm resource allocation submission file.*
```
#!/usr/bin/env bash
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --mem-per-cpu=2400
#SBATCH --job-name=thamnophis_exabayes_threading_80
#SBATCH --output=output_thamnophis_exabayes_threading_80.txt
```

Common Slurm Commands
<!-- Tables -->
| Description  | Command |
| -------------------|-------------------|
| Run analysis | sbatch [file name] |
| Pulls job ids | squeue -u [account name] |
| Check status  info on cpus | scontrol show job -d [job id] |
| Kill job | scancel [job id] |
| Check CPU load and memory  | sinfo --nodes=cpu-[number] --Format="CPUsLoad,FreeMem"|

