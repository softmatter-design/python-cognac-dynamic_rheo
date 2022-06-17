# Values used for 

# 使用するCognacのバージョン
ver_Cognac = "cognac112"
# シミュレーションに使用するコア数
core = 1

# deform_mode
deform_mode = ""
# sweep_mode
sweep_mode = ''
StrainSweepConditions = []
FrequencySweepConditions = []
temperature = 0.
min_freq = 0.
max_freq = 0.
frequency = 0.
min_strain = 0.
max_strain = 0.
strain = 0.
data_per_digit = 5

job_dir = []

# SymulationParameters
cycles = 20
# １サイクル当たりの出力回数（Step_minを割り切れる数字）
output_per_cycle = 20
# 
skipcycles = 10

# １サイクル当たりのステップ数の最小値
# (maximum value of freq*dt) = 1.0/Step_min
Step_min = int(output_per_cycle*5)
# シミュレーション時の最大 time step
dt_max = 0.01
