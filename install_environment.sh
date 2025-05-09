#!/bin/bash
mkdir -p kz_env
echo "Unpacking kz.tar.gz"
tar -xzf kz.tar.gz -C kz_env
echo "Finished unpacking kz.tar.gz"
echo "run 'source kz_env/bin/activate' in terminal to launch environment"
echo "run 'source kz_env/bin/deactivate' in terminal to leave environment"
