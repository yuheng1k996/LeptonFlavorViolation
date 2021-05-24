date
rm -rf `find -type d -name .ipynb_checkpoints`

python plots/plot_g_vs_g_0.py storage/results_g_g_0/0_max storage/results_g_g_0/0_med storage/results_g_g_0/0_min

python plots/plot_g_vs_g_1.py storage/results_g_g_1/1_max storage/results_g_g_1/1_med storage/results_g_g_1/1_min

python plots/plot_g_vs_g_2.py storage/results_g_g_2/2_max storage/results_g_g_2/2_med storage/results_g_g_2/2_min

python plots/plot_g_vs_g_3.py storage/results_g_g_3/3_max storage/results_g_g_3/3_med storage/results_g_g_3/3_min

python plots/plot_g_vs_g_4.py storage/results_g_g_4/4_max storage/results_g_g_4/4_med storage/results_g_g_4/4_min
date


