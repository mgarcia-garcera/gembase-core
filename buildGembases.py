import os
import sys
import argparse
# list of required directories



def ensure_directories(base_dir, subdirs):
    for subdir in subdirs:
        full_path = os.path.join(base_dir, subdir)
        if not os.path.exists(full_path):
            os.makedirs(full_path)
            print(f"Created directory: {full_path}")
        else:
            print(f"Directory already exists: {full_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build the gembases Directory structure")
    parser.add_argument("-d", "--dir", required=True, help="Input base directory")
    parser.add_argument("-f","--force", action="store_true", help="Force overwrite of base directory if it exists")
    parser.add_argument("-s", "--subdirectories", nargs="*", default=['Replicons','Genes','Proteins','LSTINFO','RNA'], help="List of subdirectories to create. Default is 5 standard directories.")
    
    args = parser.parse_args()
    force=args.force
    if force:
            print(f"[INFO] Removing existing base directory: {base_dir}")
            shutil.rmtree(base_dir)
            shutil.mktree(base_dir)
    ensure_directories(args.dir, args.subdirectories)
    

    