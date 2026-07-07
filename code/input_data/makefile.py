#!/usr/bin/env python3

import os
import math
import sys
import contextlib

# cal_WT_SoilMoisture.pl を Python に変換した関数
def cal_WT_SoilMoisture(point, start_year, end_year):
    """
    Calculates water temperature and soil moisture ratio based on various environmental parameters.
    """
    print("cal_WT_SoilMoisture の実行を開始します...")
    t1 = os.times()[0]
    year_length = end_year - start_year + 1
    W_soil_pre = 0
    W_snow_pre = 0

    for year in range(start_year, start_year + year_length):
        output_dir = f"./{point}/{year}/output_cal_WT"
        gnuplot_dir = f"./{point}/{year}/output_cal_WT/for_gnuplot"
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(gnuplot_dir, exist_ok=True)

        try:
            with contextlib.ExitStack() as stack:
                elv_file = stack.enter_context(open(f"./{point}/local_data/elevation_{point}.csv", 'r'))
                soilcap_file = stack.enter_context(open(f"./{point}/local_data/local_soil_capacity_{point}.csv", 'r'))
                pre_file = stack.enter_context(open(f"./{point}/{year}/input_cal_WT/{point}{year}_pre.csv", 'r'))
                ctmp_file = stack.enter_context(open(f"./{point}/{year}/input_cal_WT/{point}{year}_tmp.csv", 'r'))
                cst_file = stack.enter_context(open(f"./{point}/{year}/input_cal_WT/{point}{year}_St.csv", 'r'))
                cea_file = stack.enter_context(open(f"./{point}/{year}/input_cal_WT/{point}{year}_RH.csv", 'r'))
                cld_file = stack.enter_context(open(f"./{point}/{year}/input_cal_WT/{point}{year}_cld.csv", 'r'))
                wd_file = stack.enter_context(open(f"./{point}/{year}/input_cal_WT/{point}{year}_wind.csv", 'r'))
                ta_file = stack.enter_context(open(f"{output_dir}/{point}{year}_Ta.csv", 'w'))
                tw_file = stack.enter_context(open(f"{output_dir}/{point}{year}_Tw.csv", 'w'))
                rn_file = stack.enter_context(open(f"{output_dir}/{point}{year}_Rn.csv", 'w'))
                soilmoi_file = stack.enter_context(open(f"{output_dir}/{point}{year}_soil_moisture_ratio.csv", 'w'))
                difftwta_file = stack.enter_context(open(f"{output_dir}/{point}{year}_diffTwTa.csv", 'w'))
                roff_file = stack.enter_context(open(f"{output_dir}/{point}{year}_runoff.csv", 'w'))
                gnut_a_file = stack.enter_context(open(f"{gnuplot_dir}/{point}{year}_Ta.dat", 'w'))
                gnut_w_file = stack.enter_context(open(f"{gnuplot_dir}/{point}{year}_Tw.dat", 'w'))
                gnut_diff_file = stack.enter_context(open(f"{gnuplot_dir}/{point}{year}_Tdiff.dat", 'w'))
                gnu_rn_file = stack.enter_context(open(f"{gnuplot_dir}/{point}{year}_Rn.dat", 'w'))
                gnusoilmoi_file = stack.enter_context(open(f"{gnuplot_dir}/{point}{year}_soil_moisture_ratio.dat", 'w'))
                gnu_l_file = stack.enter_context(open(f"{gnuplot_dir}/{point}{year}_L.dat", 'w'))
                gnuroff_file = stack.enter_context(open(f"{gnuplot_dir}/{point}{year}_runoff.dat", 'w'))

                ctmp_lines = ctmp_file.readlines()
                ea_lines = cea_file.readlines()
                cst_lines = cst_file.readlines()
                pre_lines = pre_file.readlines()
                cloud_lines = cld_file.readlines()
                wind_lines = wd_file.readlines()
                elv_lines = elv_file.readlines()
                soil_lines = soilcap_file.readlines()

                for i, ctmp_line in enumerate(ctmp_lines):
                    try:
                        ctmp_data = ctmp_line.strip().split(',')
                        lon, lat = float(ctmp_data[0]), float(ctmp_data[1])
                        number_ctmp = len(ctmp_data)

                        ea_data = ea_lines[i].strip().split(',')
                        cst_data = cst_lines[i].strip().split(',')
                        pre_data = pre_lines[i].strip().split(',')
                        cloud_data = cloud_lines[i].strip().split(',')
                        wind_data = wind_lines[i].strip().split(',')
                        elv_data = elv_lines[0].strip().split(',')
                        elv = float(elv_data[2])
                        abslati = abs(lat)
                        soil_capa = soil_lines[0].strip().split(',')
                        soil_capacity = float(soil_capa[2]) * 10

                    except (ValueError, IndexError) as e:
                        print(f"警告: {point}{year} のデータファイル内に不正なデータまたは不完全な行が検出されました。この行をスキップします。エラー: {e}")
                        continue

                    ta_file.write(f"{lon},{lat},")
                    tw_file.write(f"{lon},{lat},")
                    rn_file.write(f"{lon},{lat},")
                    soilmoi_file.write(f"{lon},{lat},")
                    difftwta_file.write(f"{lon},{lat},")
                    roff_file.write(f"{lon},{lat},")

                    # === 修正箇所 ===
                    # 最初の年のみ初期値を設定し、それ以降の年は前年の値を引き継ぐ
                    if year == start_year:
                        W_soil_pre = soil_capacity * 0.8
                        W_snow_pre = 0
                    # ==================

                    for day in range(1, number_ctmp - 1):
                        try:
                            Tmp_a = float(ctmp_data[day + 1])
                            pre = float(pre_data[day + 1])
                            St0 = float(cst_data[day + 1])
                            RH = float(ea_data[day + 1])
                            cld = float(cloud_data[day + 1]) / 100.0
                            wind = float(wind_data[day + 1])
                        except (ValueError, IndexError):
                            ta_file.write("-9999,")
                            tw_file.write("-9999,")
                            rn_file.write("-9999,")
                            difftwta_file.write("-9999,")
                            soilmoi_file.write("-9999,")
                            roff_file.write("-9999,")
                            gnut_a_file.write(f"{day}\t-9999\n")
                            gnut_w_file.write(f"{day}\t-9999\n")
                            gnut_diff_file.write(f"{day}\t-9999\n")
                            gnu_rn_file.write(f"{day}\t-9999\n")
                            gnusoilmoi_file.write(f"{day}\t-9999\n")
                            gnu_l_file.write(f"{day}\t-9999\n")
                            gnuroff_file.write(f"{day}\t-9999\n")
                            # ターミナルに警告メッセージを表示
                            print(f"WARNING: {year}の{day}行目で不正なデータが検出されました。ファイルに -9999 が書き込まれました。\n")
                            continue

                        if -100 > Tmp_a or Tmp_a > 100 or pre == -9999 or St0 == -9999 or wind == -9999 or RH == -9999:
                            ta_file.write("-9999,")
                            tw_file.write("-9999,")
                            rn_file.write("-9999,")
                            difftwta_file.write("-9999,")
                            soilmoi_file.write("-9999,")
                            roff_file.write("-9999,")
                            gnut_a_file.write(f"{day}\t-9999\n")
                            gnut_w_file.write(f"{day}\t-9999\n")
                            gnut_diff_file.write(f"{day}\t-9999\n")
                            gnu_rn_file.write(f"{day}\t-9999\n")
                            gnusoilmoi_file.write(f"{day}\t-9999\n")
                            gnu_l_file.write(f"{day}\t-9999\n")
                            gnuroff_file.write(f"{day}\t-9999\n")
                            
                        # ターミナルに警告メッセージを表示
                            #print(f"WARNING: processing.pyの140行目で不正なデータが検出されました。ファイルに -9999 が書き込まれました。")
                            continue

                        if RH < 0 or RH > 100:
                            ea = (4.54 * math.exp(0.068 * Tmp_a)) / 1.332
                        else:
                            SatEvapA0 = 0.6108 * math.exp((17.27 * Tmp_a) / (Tmp_a + 237.3))
                            ea = SatEvapA0 * (RH / 100.0) * 10 / 1.332

                        SoilSurfaceTemp = Tmp_a + 0.0065 * elv
                        SoilSurfaceTemp_K = SoilSurfaceTemp + 273.15

                        if SoilSurfaceTemp_K * elv == 0:
                            pressure_modi = 1013.25
                        else:
                            pressure_modi = 1013.25 * ((1 - 0.0065 / SoilSurfaceTemp_K * elv) ** 5.2552)

                        fdf = -0.1742 + 0.00116 * pressure_modi
                        if fdf <= 0:
                            fdf = 1.0

                        if W_soil_pre / soil_capacity * 100 >= 30 and Tmp_a > 5:
                            LAI_data = 2.0
                        else:
                            LAI_data = 0.01

                        Fn_from_Tw = (Tmp_a + 273.15) ** 4 * 0.0000000567
                        ea_fdf = ea * 1.332 * fdf
                        Fnair_coff = 46.5 * (ea_fdf / (Tmp_a + 273.15))
                        datasq = 1.2 + 3 * Fnair_coff
                        Fnair_CearSky_C = 1 - (1 + Fnair_coff) * math.exp(-math.sqrt(datasq)) if datasq > 0 else 0
                        Fnair_CearSky = Fnair_CearSky_C * 0.0000000567 * (Tmp_a + 273.15) ** 4
                        Fnair_coff_C = 1 - (0.095 - 0.0006 * ea_fdf) * cld - (0.66 - 0.0044 * ea_fdf) * cld
                        Fn_air = 0.0000000567 * (Tmp_a + 273.15) ** 4 * (1 - (1 - (Fnair_CearSky / (0.0000000567 * (Tmp_a + 273.15) ** 4))) * Fnair_coff_C)

                        albedo_warter = 0.06
                        if Tmp_a <= 0:
                            albedo_warter = 0.70

                        if LAI_data > 0 and Tmp_a > 0:
                            alb_plant = 0.22 - (0.22 - 0.06) * math.exp(-0.5 * LAI_data)
                            St_abave_plant = St0 * math.exp(-0.5 * LAI_data)
                            Stp = St_abave_plant * (1 - alb_plant)
                        else:
                            Stp = St0 * (1 - albedo_warter)

                        Fn_from_plant = (Tmp_a + 273.15) ** 4 * 0.0000000567 * (1 - math.exp(-0.5 * LAI_data)) / 2 if LAI_data > 0 else 0
                        Rn = Stp + Fn_air + Fn_from_plant - Fn_from_Tw

                        h_cal = 0.0002
                        h_k = 2
                        Rn3 = Rn * 2.06 / 60 / 60 / 24
                        hR = 4 * 0.00000000000136 * ((Tmp_a + 273.15) ** 3)
                        CO = h_k * (4.75 - ea) - (Tmp_a) * (1 + (hR / h_cal)) - (Rn3 / h_cal)
                        C1 = h_k * 0.265 + 1 + (hR / h_cal)
                        C2 = 0.01876
                        C3 = 0.0008032
                        a = C2 / C3
                        b = C1 / C3
                        c = CO / C3
                        p = b / 3 - a ** 2 / 9
                        q = c - (a * b) / 3 + (2 * (a ** 3)) / 27
                        q24p = q ** 2 + 4 * p ** 3
                        if q24p < 0:
                            q24p = 0
                        alfa = (1 / 2) * (-q + math.sqrt(q24p))
                        t_val = -math.sqrt(q24p)
                        beta = -((1 / 2) * (-q + t_val))
                        estimatedWT = alfa ** (1 / 3) - beta ** (1 / 3) - a / 3

                        Tmp_w = estimatedWT
                        diff_TwTa = Tmp_w - Tmp_a

                        ta_file.write(f"{Tmp_a},")
                        tw_file.write(f"{Tmp_w},")
                        rn_file.write(f"{Rn:.2f},")
                        difftwta_file.write(f"{diff_TwTa},")

                        gnut_a_file.write(f"{day}\t{Tmp_a}\n")
                        gnut_w_file.write(f"{day}\t{Tmp_w}\n")
                        gnut_diff_file.write(f"{day}\t{diff_TwTa}\n")
                        gnu_rn_file.write(f"{day}\t{Rn:.2f}\n")
                        gnu_l_file.write(f"{day}\t{Fn_air}\n")

                        diff_a = 4.75
                        diff_b = 0.265
                        diff_c = 0.00938
                        diff_d = 0.0004016
                        es = (diff_a + diff_b * Tmp_a + diff_c * (Tmp_a ** 2) + diff_d * (Tmp_a ** 3))
                        e_diff = es - ea
                        e_diff_kPa = e_diff / 7.50062
                        fix_e_diff_kPa = e_diff_kPa * 1.0
                        if fix_e_diff_kPa < 0:
                            fix_e_diff_kPa = 0

                        min_tmp, max_tmp = Tmp_a - 0.01, Tmp_a + 0.01
                        min_es = (diff_a + diff_b * min_tmp + diff_c * (min_tmp ** 2) + diff_d * (min_tmp ** 3))
                        max_es = (diff_a + diff_b * max_tmp + diff_c * (max_tmp ** 2) + diff_d * (max_tmp ** 3))
                        delta = (max_es / 7.50062 - min_es / 7.50062) / (max_tmp - min_tmp)
                        Pa = pressure_modi / 10
                        gamma = (6.62 * (10 ** (-4))) * Pa * 1.006
                        fix_Rn = Rn * (60 * 60 * 24) / (10 ** 6)
                        G = 0
                        U = wind * 0.8
                        if U < 0:
                            U = 2
                        penman_bunshi = 0.408 * delta * (fix_Rn - G) + gamma * (900 / (Tmp_a + 273)) * U * fix_e_diff_kPa
                        penman_bunbo = delta + gamma * (1 + 0.34 * U)
                        ETo = penman_bunshi / penman_bunbo
                        if ETo < 0:
                            ETo = 0

                        P_rain = 0
                        P_snow = 0
                        melt = 0
                        runoff = 0

                        if Tmp_a < -1:
                            P_snow = pre
                            P_rain = 0
                            melt_val = 2.63 + 2.55 * Tmp_a + 0.0912 * Tmp_a * P_rain
                            melt = max(0, min(melt_val, W_snow_pre + P_snow))
                            W_snow = W_snow_pre + P_snow - melt
                        else:
                            P_snow = 0
                            P_rain = pre
                            if W_snow_pre > 0:
                                melt_val = 2.63 + 2.55 * Tmp_a + 0.0912 * Tmp_a * P_rain
                                melt = max(0, min(melt_val, W_snow_pre + P_snow))
                            else:
                                melt = 0
                            W_snow = W_snow_pre + P_snow - melt

                        depth = 0.8
                        Sa = 200
                        dep = 0.5
                        soil_a = (W_soil_pre + P_rain + melt) * depth
                        soil_b = Sa * depth * (1 - dep)

                        if soil_a >= soil_b:
                            ETa = ETo
                        else:
                            rho = (W_soil_pre + P_rain + melt) / (Sa * (1 - dep))
                            ETa = ETo * rho

                        W_soil_1 = W_soil_pre + P_rain + melt - ETa
                        W_soil_2 = soil_capacity

                        if W_soil_1 < W_soil_2:
                            W_soil = W_soil_1
                        else:
                            W_soil = W_soil_2

                        runoff = W_soil_1 - W_soil
                        roff_file.write(f"{runoff},")
                        gnuroff_file.write(f"{day}\t{runoff}\n")

                        W_snow_pre = W_snow
                        W_soil_pre = W_soil
                        soil_moisture = (W_soil_pre / soil_capacity) * 100
                        soilmoi_file.write(f"{soil_moisture},")
                        gnusoilmoi_file.write(f"{day}\t{soil_moisture}\n")

                ta_file.write("\n")
                tw_file.write("\n")
                rn_file.write("\n")
                soilmoi_file.write("\n")
                difftwta_file.write("\n")
                roff_file.write("\n")

        except FileNotFoundError as e:
            print(f"エラー: {e}. cal_WT_SoilMoisture の実行に必要な入力ファイルがすべて存在することを確認してください。")
            sys.exit(1)

    t2 = os.times()[0]
    t3 = t2 - t1
    print(f"\ncal_WT_SoilMoisture の実行が完了しました。所要時間: {t3} 秒\n")

# make_Dchange.pl を Python に変換した関数
def make_Dchange(validation_point, start_year, end_year):
    """
    Calculates the daily change in day length and saves the data to a CSV file.
    """
    print("make_Dchange の実行を開始します...")
    
    # 複数年のデータを結合して処理
    all_photp = []
    year_data_counts = {}
    
    for year in range(start_year, end_year + 1):
        try:
            with open(f"./{validation_point}/{validation_point}_photp.csv", 'r') as f:
                line = f.readlines()[0].strip().split(',')
                photp_data = [float(d) for d in line[2:] if d.strip()]
                all_photp.extend(photp_data)
                year_data_counts[year] = len(photp_data)
        except (FileNotFoundError, ValueError, IndexError) as e:
            print(f"致命的なエラー: {validation_point}_photp.csv の読み込み中にエラーが発生しました: {e}")
            sys.exit(1)

    total_data_index = 0
    for year in range(start_year, end_year + 1):
        output_dir = f"./{validation_point}/{year}"
        os.makedirs(output_dir, exist_ok=True)
        
        with open(f"{output_dir}/{validation_point}{year}_Dch.csv", 'w') as dch_file:
            dch_file.write("0,1")
            
            for day_counter in range(1, year_data_counts[year]):
                dch_file.write(f",{all_photp[total_data_index + day_counter] - all_photp[total_data_index + day_counter - 1]}")
            
            dch_file.write("\n")
        
        total_data_index += year_data_counts[year]
    
    print("\nmake_Dchange の実行が完了しました。\n")

# make_Dweek.pl を Python に変換した関数
def make_Dweek(validation_point, start_year, end_year):
    """
    Calculates the 7-day average of day length and saves the data to a CSV file.
    """
    print("make_Dweek の実行を開始します...")
    
    all_photp = []
    year_data_counts = {}
    
    for year in range(start_year, end_year + 1):
        try:
            with open(f"./{validation_point}/{validation_point}_photp.csv", 'r') as f:
                line = f.readlines()[0].strip().split(',')
                photp_data = [float(d) for d in line[2:] if d.strip()]
                all_photp.extend(photp_data)
                year_data_counts[year] = len(photp_data)
        except (FileNotFoundError, ValueError, IndexError) as e:
            print(f"致命的なエラー: {validation_point}_photp.csv の読み込み中にエラーが発生しました: {e}")
            sys.exit(1)
    
    total_data_index = 0
    for year in range(start_year, end_year + 1):
        output_dir = f"./{validation_point}/{year}"
        os.makedirs(output_dir, exist_ok=True)
        
        with open(f"{output_dir}/{validation_point}{year}_Dweek.csv", 'w') as dwk_file:
            dwk_file.write("0,1")
            
            for day_counter in range(year_data_counts[year]):
                photp_sum = 0
                count = 0
                for d in range(7):
                    if total_data_index + day_counter - d >= 0:
                        photp_sum += all_photp[total_data_index + day_counter - d]
                        count += 1
                
                dweek = photp_sum / max(1, count)
                dwk_file.write(f",{dweek}")
            
            dwk_file.write("\n")
        
        total_data_index += year_data_counts[year]
    
    print("\nmake_Dweek の実行が完了しました。\n")


# make_Taweek.pl を Python に変換した関数
def make_Taweek(validation_point, start_year, end_year):
    """
    Calculates the 7-day average of air temperature and saves the data to a CSV file.
    """
    print("make_Taweek の実行を開始します...")
    
    all_ta_data = []
    year_data_counts = {}
    
    for year in range(start_year, end_year + 1):
        try:
            with open(f"./{validation_point}/{year}/output_cal_WT/{validation_point}{year}_Ta.csv", 'r') as f:
                line = f.readlines()[0].strip().split(',')
                ta_data = [float(d) if d.strip() and d.strip() != '-9999' else -9999 for d in line[2:]]
                all_ta_data.extend(ta_data)
                year_data_counts[year] = len(ta_data)
        except (FileNotFoundError, ValueError, IndexError) as e:
            print(f"致命的なエラー: {validation_point}{year}_Ta.csv の読み込み中にエラーが発生しました: {e}")
            sys.exit(1)

    total_data_index = 0
    for year in range(start_year, end_year + 1):
        output_dir = f"./{validation_point}/{year}"
        os.makedirs(output_dir, exist_ok=True)
        
        with open(f"{output_dir}/{validation_point}{year}_Taweek.csv", 'w') as taw_file:
            taw_file.write("0,1")
            
            for day_counter in range(year_data_counts[year]):
                ta_sum = 0
                count = 0
                for d in range(7):
                    if total_data_index + day_counter - d >= 0 and all_ta_data[total_data_index + day_counter - d] != -9999:
                        ta_sum += all_ta_data[total_data_index + day_counter - d]
                        count += 1
                
                taweek = ta_sum / max(1, count) if count > 0 else -9999
                taw_file.write(f",{taweek}")
            
            taw_file.write("\n")
        
        total_data_index += year_data_counts[year]
    
    print("\nmake_Taweek の実行が完了しました。\n")

# make_Twweek.pl を Python に変換した関数
def make_Twweek(validation_point, start_year, end_year):
    """
    Calculates the 7-day average of water temperature and saves the data to a CSV file.
    """
    print("make_Twweek の実行を開始します...")
    
    all_tw_data = []
    year_data_counts = {}
    
    for year in range(start_year, end_year + 1):
        try:
            with open(f"./{validation_point}/{year}/output_cal_WT/{validation_point}{year}_Tw.csv", 'r') as f:
                line = f.readlines()[0].strip().split(',')
                tw_data = [float(d) if d.strip() and d.strip() != '-9999' else -9999 for d in line[2:]]
                all_tw_data.extend(tw_data)
                year_data_counts[year] = len(tw_data)
        except (FileNotFoundError, ValueError, IndexError) as e:
            print(f"致命的なエラー: {validation_point}{year}_Tw.csv の読み込み中にエラーが発生しました: {e}")
            sys.exit(1)
    
    total_data_index = 0
    for year in range(start_year, end_year + 1):
        output_dir = f"./{validation_point}/{year}"
        os.makedirs(output_dir, exist_ok=True)
        
        with open(f"{output_dir}/{validation_point}{year}_Twweek.csv", 'w') as tww_file:
            tww_file.write("0,1")
            
            for day_counter in range(year_data_counts[year]):
                tw_sum = 0
                count = 0
                for d in range(7):
                    if total_data_index + day_counter - d >= 0 and all_tw_data[total_data_index + day_counter - d] != -9999:
                        tw_sum += all_tw_data[total_data_index + day_counter - d]
                        count += 1
                
                twweek = tw_sum / max(1, count) if count > 0 else -9999
                tww_file.write(f",{twweek}")
            
            tww_file.write("\n")
        
        total_data_index += year_data_counts[year]
    
    print("\nmake_Twweek の実行が完了しました。\n")

# make_data_for_simulation2.pl を Python に変換した関数
def make_data_for_simulation2(validation_point, start_year, end_year):
    """
    Combines climate data from various CSV files into a single output file.
    """
    print("make_data_for_simulation2 の実行を開始します...")
    total_count = 0

    for year in range(start_year, end_year + 1):
        try:
            output_dir = f"./{validation_point}/{year}"
            os.makedirs(output_dir, exist_ok=True)

            with contextlib.ExitStack() as stack:
                cclim_file = stack.enter_context(open(f"{output_dir}/{validation_point}{year}_climdata2.csv", 'w'))
                at_file = stack.enter_context(open(f"{output_dir}/input_cal_WT/{validation_point}{year}_tmp.csv", 'r'))
                wt_file = stack.enter_context(open(f"{output_dir}/output_cal_WT/{validation_point}{year}_Tw.csv", 'r'))
                dwk_file = stack.enter_context(open(f"{output_dir}/{validation_point}{year}_Dweek.csv", 'r'))
                sm_file = stack.enter_context(open(f"{output_dir}/output_cal_WT/{validation_point}{year}_soil_moisture_ratio.csv", 'r'))
                pre_file = stack.enter_context(open(f"{output_dir}/input_cal_WT/{validation_point}{year}_pre.csv", 'r'))
                roff_file = stack.enter_context(open(f"{output_dir}/output_cal_WT/{validation_point}{year}_runoff.csv", 'r'))
                dch_file = stack.enter_context(open(f"{output_dir}/{validation_point}{year}_Dch.csv", 'r'))
                taw_file = stack.enter_context(open(f"{output_dir}/{validation_point}{year}_Taweek.csv", 'r'))
                tww_file = stack.enter_context(open(f"{output_dir}/{validation_point}{year}_Twweek.csv", 'r'))

                at_lines = at_file.readlines()
                wt_lines = wt_file.readlines()
                dwk_lines = dwk_file.readlines()
                sm_lines = sm_file.readlines()
                pre_lines = pre_file.readlines()
                roff_lines = roff_file.readlines()
                dch_lines = dch_file.readlines()
                taw_lines = taw_file.readlines()
                tww_lines = tww_file.readlines()

                # すべてのファイルで最も短い行の長さを取得し、それをもとにループを回す
                file_lines = [at_lines, wt_lines, dwk_lines, sm_lines, pre_lines, roff_lines, dch_lines, taw_lines, tww_lines]
                min_lines = min(len(lines) for lines in file_lines)

                for i in range(min_lines):
                    # 各行を読み込み、データの数を比較して最小値に合わせる
                    data_at = at_lines[i].strip().split(',')
                    data_wt = wt_lines[i].strip().split(',')
                    data_dwk = dwk_lines[i].strip().split(',')
                    data_sm = sm_lines[i].strip().split(',')
                    data_pre = pre_lines[i].strip().split(',')
                    data_roff = roff_lines[i].strip().split(',')
                    data_dch = dch_lines[i].strip().split(',')
                    data_taw = taw_lines[i].strip().split(',')
                    data_tww = tww_lines[i].strip().split(',')

                    # 各ファイルの日ごとのデータの長さを取得
                    data_lists = [data_at, data_wt, data_dwk, data_sm, data_pre, data_roff, data_dch, data_taw, data_tww]
                    min_data_len = min(len(d) for d in data_lists)

                    # 最小のデータ長に合わせてループを回す
                    for j in range(2, min_data_len):
                        try:
                            # 結合して出力ファイルに書き込む
                            cclim_file.write(f"{data_at[j]},{data_wt[j]},{data_dwk[j]},{data_sm[j]},{data_pre[j]},{data_roff[j]}")
                            if j < min_data_len - 1:
                                cclim_file.write("\n")
                            total_count += 1
                        except IndexError as e:
                            print(f"警告: {validation_point}{year} のデータ結合中に不整合なデータが検出されました。この行をスキップします。エラー: {e}")
                            continue

                cclim_file.write("\n")

        except FileNotFoundError as e:
            print(f"エラー: {e}. make_data_for_simulation2 の実行に必要な入力ファイルがすべて存在することを確認してください。")
            sys.exit(1)

    print(f"\nmake_data_for_simulation2 の実行が完了しました。\n")
    print(f"合計 {total_count} 件のデータを結合しました。")


# プログラムのメインエントリポイント
if __name__ == "__main__":
    try:
        # ユーザー入力
        print("--- プログラム実行のための情報を入力してください ---")
        point = input("処理する都市名を入力してください (例: tokyo): ")
        start_year_str = input("開始年を入力してください (例: 2002): ")
        end_year_str = input("終了年を入力してください (例: 2013): ")

        # 文字列を整数に変換
        start_year = int(start_year_str)
        end_year = int(end_year_str)
        year_length = end_year - start_year + 1

        if start_year > end_year:
            print("エラー: 開始年は終了年より後の年を指定することはできません。")
            sys.exit(1)

        print(f"\n--- {point} の {start_year} 年から {end_year} 年までの処理を開始します ---\n")

        # 実行順序: 各ステップは前のステップの出力に依存しています。
        cal_WT_SoilMoisture(point, start_year, end_year)
        make_Dchange(point, start_year, end_year)
        make_Dweek(point, start_year, end_year)
        make_Taweek(point, start_year, end_year)
        make_Twweek(point, start_year, end_year)
        make_data_for_simulation2(point, start_year, end_year)

        print("\n--- 全プログラムの実行が完了しました ---\n")

    except ValueError:
        print("エラー: 年は半角数字で入力してください。")
        sys.exit(1)
    except Exception as e:
        print(f"致命的なエラーが発生しました: {e}")
        sys.exit(1)