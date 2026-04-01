import os
import requests
from tqdm import tqdm
import argparse

def download_isform(symbols, data_dir):
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        print(f"Created directory: {data_dir}")
    base_url = "https://zenodo.org/records/19365039/files"
    print("--- Downloading Files ---")
    for symbol in symbols:
            if symbol == "T.aestivum":
                files_to_download = [
                    f"{symbol}_isoform.A.fa", 
                    f"{symbol}_isoform.B.fa", 
                    f"{symbol}_isoform.D.fa"
                ]
            else:
                files_to_download = [f"{symbol}_isoform.fa"]
            for filename in files_to_download:
                target_path = os.path.join(data_dir, filename)
                url = f"{base_url}/{filename}?download=1"
                if os.path.exists(target_path):
                    print(f"File {filename} already exists, skipping...")
                    continue
                print(f"Downloading {filename}...")
                try:
                    with requests.get(url, stream=True, timeout=60) as r:
                        r.raise_for_status()
                        with open(target_path, 'wb') as f:
                            for chunk in r.iter_content(chunk_size=8192):
                                f.write(chunk)
                    print(f"Successfully downloaded {filename}")
                except Exception as e:
                    print(f"Error downloading {filename}: {e}")
                    return False
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CellBlaster Data Downloader")
    parser.add_argument("-s", "--symbols", nargs='+', 
                        default=["A.thaliana", "G.max", "L.japonicus"],
                        help="List of species symbols to download (e.g., -s A.thaliana G.max)")
    parser.add_argument("-o", "--data_dir", default="./Download_FASTA",
                        help="Directory to save the downloaded files")
    args = parser.parse_args()
    symbols = args.symbols
    data_dir = os.path.abspath(os.path.expanduser(args.data_dir))
    os.makedirs(data_dir, exist_ok=True)
    os.chdir(data_dir)
    # --- ownloading Files ---
    print("--- Step 1: Downloading Files ---")
    download_success = download_isform(symbols, data_dir)