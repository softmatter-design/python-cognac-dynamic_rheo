# Values

# 使用するCognacのバージョン
ver_Cognac = "cognac112"
# シミュレーションに使用するコア数
core = 1

read_udf = ''
deform_mode = ''
sweep_mode = ''
Conditions = []
base_dir = ''
base_udf = 'base_out.udf'
Sim_conditions = []
target_dir = ''
subdir_list = []

# SymulationParameters
total_cycles = 25
# １サイクル当たりの出力回数（Step_minを割り切れる数字）
output_cycle = 20
# 
skipcycles = 5
# １サイクル当たりのステップ数の最小値 
step_min = int(output_cycle*5)
# シミュレーション時の最大 time step
dt_max = 0.01
time = []

batch = ''
restart_udf = ''