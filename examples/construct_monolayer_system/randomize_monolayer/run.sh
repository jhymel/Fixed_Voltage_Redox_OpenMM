#!/bin/bash

job=$(qsub simulated_heating.pbs)
qsub -W depend=afterok:${job} simulated_annealing.pbs

