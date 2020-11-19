#!/bin/bash
srun -c 1 --mem 100G -p gpu --gres=gpu:1 --pty bash -i
