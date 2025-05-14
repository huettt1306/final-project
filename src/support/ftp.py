import json
from ftplib import FTP

ftp_server = 'ftp.sra.ebi.ac.uk'
ftp = FTP(ftp_server)
ftp.login()

base_path = '/vol1/run/'
result_path = '/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/test.txt'

with open("/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/trio.json", "r") as file:
        sample_data = json.load(file)

sample_names = set()
for sample_group in sample_data.values():
    sample_names.update(sample_group.values())  

print (sample_names)

found_samples = set()

target_dirs = {
    'ERR324': (3242356, 3242356),  
}

def check_and_save_paths():
    with open(result_path, 'w') as f_out:
        for dir_a, (start_number, end_number) in target_dirs.items():
            print(f"Đang duyệt thư mục chính: {dir_a}")
            ftp.cwd(base_path)
            ftp.cwd(dir_a)

            for current_number in range(start_number, end_number + 1):
                dir_b = f"ERR{current_number:07d}" 
                print(dir_b)
                try:
                    ftp.cwd(dir_b) 
                    for filename in ftp.nlst():
                        if filename.endswith(".final.cram"):
                            parts = filename.split('.')
                            if parts:
                                sample_name = parts[0]  
                                print(sample_name)
                                if sample_name in sample_names:
                                    remote_path = f"{ftp_server}{base_path}{dir_a}/{dir_b}/{filename}"
                                    f_out.write(f"{remote_path}\n")
                                    found_samples.add(sample_name)  
                            else:
                                print(f"Bỏ qua file không hợp lệ: {filename}")
                    ftp.cwd(base_path + dir_a)  
                except Exception as e:
                    print(f"Không tìm thấy thư mục {dir_b}, dừng duyệt.")
                    break

                current_number += 1  

    print("Đã lưu các đường dẫn vào ftppath.txt.")


check_and_save_paths()

missing_samples = sample_names - found_samples
if missing_samples:
    print("Các mẫu chưa có đường dẫn tải:")
    for missing in missing_samples:
        print(missing)
else:
    print("Tất cả các mẫu đều đã có đường dẫn tải.")

ftp.quit()
