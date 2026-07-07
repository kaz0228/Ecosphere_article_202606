import pandas as pd
import os
import numpy as np

def process_data(city, year):
    """
    指定された都市と年のCSVファイルから複数の列のデータを抽出し、
    新しい形式でそれぞれ別のファイルに保存する関数。

    Args:
        city (str): 都市名。
        year (int): 年。
    """
    input_file_path = f'./{city}/row/{city}_{year}.csv'
    output_dir_path = f'./{city}/processing/{year}/'

    # 入力ファイルが存在するか確認
    if not os.path.exists(input_file_path):
        print(f"エラー: 指定されたファイルが見つかりません - {input_file_path}")
        return

    # 出力ディレクトリが存在しない場合は作成
    os.makedirs(output_dir_path, exist_ok=True)

    try:
        # CSVファイルを読み込み、必要なすべての列を抽出
        # 列のインデックスは0から始まるため、2, 5, 9, 12, 15, 18列目はそれぞれ1, 4, 8, 11, 14, 17となる
        df = pd.read_csv(
            input_file_path,
            header=None,
            skiprows=5,
            usecols=[1, 4, 8, 11, 14, 17],
            encoding='cp932'
        )

        # 抽出した各列を個別のDataFrameに分割
        df_tmp = df.iloc[:, 0].to_frame()
        df_pre = df.iloc[:, 1].to_frame()
        df_rh = df.iloc[:, 2].to_frame()
        df_st = df.iloc[:, 3].to_frame()
        df_cld = df.iloc[:, 4].to_frame()
        df_wind = df.iloc[:, 5].to_frame()

        # cld（雲量）の値を1/10にし、小数点以下第2位で切り捨て
        df_cld = df_cld / 10
        df_cld = np.floor(df_cld * 100) / 100
        
        # St（全天日射量）の単位をMJ/㎡からWに変換し、小数点以下第2位で四捨五入
        # 1日あたりの平均値として計算（1日 = 24時間 * 3600秒 = 86400秒）
        df_st = (df_st * 10**6) / 86400
        df_st = df_st.round(2)

        # 各DataFrameに0と1の行を挿入し、転置して保存する汎用関数
        def save_processed_data(data_frame, file_name):
            new_rows = pd.DataFrame({data_frame.columns[0]: [0, 1]})
            processed_df = pd.concat([new_rows, data_frame], ignore_index=True)
            processed_df_transposed = processed_df.T
            output_file_path = f'{output_dir_path}{file_name}'
            processed_df_transposed.to_csv(output_file_path, index=False, header=False, encoding='cp932')
            print(f"データが正常に処理され、'{output_file_path}'に保存されました。")

        # 各ファイルを保存
        save_processed_data(df_tmp, f'{city}{year}_tmp.csv')
        save_processed_data(df_pre, f'{city}{year}_pre.csv')
        save_processed_data(df_rh, f'{city}{year}_RH.csv')
        save_processed_data(df_st, f'{city}{year}_St.csv')
        save_processed_data(df_cld, f'{city}{year}_cld.csv')
        save_processed_data(df_wind, f'{city}{year}_wind.csv')
            
    except UnicodeDecodeError:
        print("エンコーディングエラー: 'cp932'でもファイルをデコードできませんでした。ファイルのエンコーディングを確認してください。")
    except Exception as e:
        print(f"処理中に予期せぬエラーが発生しました: {e}")

# プログラムの開始地点
if __name__ == "__main__":
    # ユーザーからの入力を受け取る
    city_input = input("処理する都市名を入力してください (例: tokyo): ")
    year_input = input("処理する年を入力してください (例: 2014): ")
    
    try:
        year_input = int(year_input)
        process_data(city_input, year_input)
    except ValueError:
        print("エラー: 年は数字で入力してください。")