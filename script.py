src_file1 = 'control.R'
src_file2 = 'temp['
src1 = open(src_file1, 'r', encoding='utf-8')
src2 = open(src_file2, 'r', encoding='utf-8')
outfile = open('outfile', 'w')

val_dict = {}
for line in src1.readlines():
    pos = line.find('#')
    note, key, val = '', '', ''
    if pos != -1:
        note = line[pos:]
    key, sep, val = line[0:pos].split()
    val_dict[key] = [val, note.replace('\n', '')]

for line in src2.readlines():
    key, sep, val = line.split()
    curr = val_dict.get(key, ['', ''])
    curr[0] = val
    val_dict[key] = curr

for key in val_dict:
    s = key + ' <- ' + val_dict[key][0]
    outfile.write(s.ljust(40, ' ') + val_dict[key][1] + '\n')
