import pandas as pd
import numpy as np
import os

def calculate_weekly_averages(file_path):
    """
    指定されたCSVファイルからデータを読み込み、7日間（週）ごとの平均を計算する関数。
    
    Args:
        file_path (str): データが保存されているCSVファイルへのパス。

    Returns:
        pd.DataFrame: 7日間平均が計算された新しいデータフレーム。
    """
    try:
        # ヘッダーを5行スキップしてデータを読み込み、列名を明示的に指定
        df = pd.read_csv(file_path, header=None, skiprows=5, usecols=[0, 1, 4, 8, 11, 14, 17], encoding='utf-8')

        # 列名を新しいファイル形式に合わせて変更
        df.columns = ['Date', '気温(℃)', '降水量(mm)', '日射量(MJ/㎡)', '風速(m/s)', '湿度(％)', '雲量(10分比)']
        
        # 欠損値（'--'や空白など）をNaNに変換し、数値型に変換
        df = df.replace(['--','///',' '], np.nan)
        for col in df.columns[1:]:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # 日付の列を直接変換
        df['Date'] = pd.to_datetime(df['Date'], format="%Y/%m/%d")
        df.set_index('Date', inplace=True)
        
        # 7日間の移動平均を計算し、週ごとの平均をリサンプリングで計算
        weekly_averages = df.resample('W-MON').mean()
        
        return weekly_averages

    except FileNotFoundError:
        print(f"エラー: 指定されたファイル '{file_path}' が見つかりません。")
        return None
    except Exception as e:
        print(f"データの処理中にエラーが発生しました: {e}")
        return None

if __name__ == '__main__':
    # ユーザーに入力を促す
    point = input("処理する地点名を入力してください（例：kumagaya）: ")
    year = input("処理する年を入力してください（例：2016）: ")
    
    # 修正箇所: 入力ファイルパスの変更
    input_file = f'./{point}/row/{point}_{year}.csv'
    output_dir = './week_average'
    output_file = f'{output_dir}/{point}_{year}_weekly_averages.csv'

    # 週ごとの平均を計算
    weekly_data = calculate_weekly_averages(input_file)
    
    if weekly_data is not None:
        # 出力ディレクトリが存在しない場合は作成
        os.makedirs(output_dir, exist_ok=True)

        print("週ごとの7日間平均データ:")
        print(weekly_data)
        
        # 結果を新しいCSVファイルとして出力
        weekly_data.to_csv(output_file)
        print(f"\n結果は '{output_file}' に保存されました。")