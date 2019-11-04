# 为了防止数据过小出现e-100次方导致matlab无法读取文件数据
import os
import re
file_path = '../orbit_results'
re_pattern = '-?\d\.\d*-\d*'
re_replace = '0.0000E+00'
for file_name in os.listdir(file_path):
    edit_file = ''.join(['../orbit_results/', file_name])
    with open(edit_file, 'r+', encoding='utf-8') as f:
        file_lines = f.readlines()
        i = 0
        for each_line in file_lines:
            new_line = re.sub(re_pattern, re_replace, each_line)
            file_lines[i] = new_line
            i = i + 1
    with open(edit_file, 'w', encoding='utf-8') as f:
        f.writelines(file_lines)
