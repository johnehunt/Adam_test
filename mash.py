import os
import subprocess

def run_mash(query_file, reference_file):
    # Define Mash command
    mash_cmd = ["mash", "dist", query_file, reference_file]

    # Run Mash and capture the output
    result = subprocess.run(mash_cmd, capture_output=True, text=True)

    # Parse the Mash output to get the distance
    distance = float(result.stdout.strip().split()[2])

    return distance

def main():
    # Set the path to the folder containing GBK files
    gbk_folder = "/path/to/gbk/folder"

    # Set the path to the reference biosynthetic gene cluster
    reference_bgc = "/path/to/reference.bgc"

    # Iterate through all GBK files in the folder
    for gbk_file in os.listdir(gbk_folder):
        if gbk_file.endswith(".gbk"):
            # Full path to the current GBK file
            gbk_path = os.path.join(gbk_folder, gbk_file)

            # Run Mash against the reference biosynthetic gene cluster
            distance = run_mash(gbk_path, reference_bgc)

            # Print the result
            print(f"{gbk_file}: Mash distance = {distance}")

if __name__ == "__main__":
    main()
