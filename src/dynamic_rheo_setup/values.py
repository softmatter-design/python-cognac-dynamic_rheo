# Values used for 

# 使用するCognacのバージョン
ver_Cognac = "cognac112"
# シミュレーションに使用するコア数
core = 1
# 計算条件
Method = ""
# 基準周波数
Freq_Base = 0.
# 基準ひずみ
Str_Base = 0.
# スイープレンジ
Range = [1, 0.1]
#
# 解像度
# Resol = [1.0, 0.79, 0.63, 0.5, 0.4, 0.32, 0.25, 0.2, 0.16, 0.13]
# Resol = [1.0, 0.63, 0.4, 0.25, 0.16]
Resol = [1., 0.56, 0.32, 0.18]
# Resol = [1., 0.5, 0.2]
#
# シミュレーション時の最大 time step
dt_max = 0.01
################################################################################
#========== parameters ==========
# １サイクル当たりのステップ数の最小値
# (maximum value of freq*dt) = 1.0/Step_min
Step_min = 200
# １サイクル当たりの出力回数（Step_minを割り切れる数字）
output_per_cycle = 40
# シミュレーションを行うサイクル数
num_cycles = 30
# name of analysis script
py_fname = "Analysis_V2.py"
