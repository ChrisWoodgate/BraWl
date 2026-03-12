import os
import shutil

# Function to create directories and copy only .inp files from the source directory
def create_directories_and_copy_files(brawl_type, source_dir, num_walkers, overlaps):
    # Cores list for which directories will be created
    windows = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    # Check if source directory exists
    if not os.path.exists(source_dir):
        print(f"Source directory '{source_dir}' does not exist.")
        return

    # Loop through the range of values for windows
    for brawl in brawl_type:
      for window in windows:
        for walkers in num_walkers:
          for overlap in overlaps:
            for i in range(5):
                # Create the directory name following the naming convention
                dir_name = f'{brawl}_{walkers:02d}_{window:02d}_{overlap:02d}_{i+1}'

                # Create the directory if it does not exist
                if not os.path.exists(dir_name):
                    os.makedirs(dir_name)

                for item in os.listdir(source_dir):
                    source_item = os.path.join(source_dir, item)
                    destination_item = os.path.join(dir_name, item)
                    
                    # Only copy .inp or .vij files, and skip directories
                    if os.path.isfile(source_item) and item.endswith(('.inp', '.vij')):
                        shutil.copy2(source_item, destination_item)  # copy2 preserves metadata
                
                file_path = os.path.join(dir_name, "wl_input.inp")
                with open(file_path, "r") as file:
                    lines = file.readlines()

                for j in range(len(lines)):
                  if lines[j].strip().startswith("num_windows="):
                      lines[j] = f"num_windows={window}\n"
                  elif lines[j].strip().startswith("performance="):
                      lines[j] = f"performance={brawl}\n"
                  elif lines[j].strip().startswith("bin_overlap="):
                      lines[j] = f"bin_overlap={overlap/100}\n"

                with open(file_path, "w") as file:
                    file.writelines(lines)

# Example call to the function with brawl_type and source_dir
brawl_type = [0, 1, 2, 3, 4, 5]
source_dir = '.'
num_walkers = [1]
overlaps = [0, 10, 25, 50, 75]

brawl_type = [0, 2, 4]
num_walkers = [ 2, 3, 4, 5, 6]
overlaps = [25]

create_directories_and_copy_files(brawl_type, source_dir, num_walkers, overlaps)
