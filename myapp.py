import os
import re
import sys
import subprocess
from pathlib import Path
from flask import Flask, render_template, request
from markupsafe import Markup
from werkzeug.utils import secure_filename

app = Flask(__name__)

# 将zip函数添加到jinja2的全局命名空间
app.jinja_env.globals.update(zip=zip)


app.config['UPLOAD_FOLDER'] = './uploads'


@app.route('/guidelines')
def guidelines():
    return render_template('guidelines.html')


# @app.route('/')
# def ss():
#     return render_template('result.html')

@app.route('/', methods=['GET', 'POST'])
def sirna_designer():
    if request.method == 'POST':
        # 这是一个POST请求，处理表单数据
        nom = request.form.get('nom')
        sequence = request.form.get('sequence')

        taille_motif = request.form.get('taille_motif')
        consGC = request.form.get('consGC')
        consAT = request.form.get('consAT')
        consGT = request.form.get('consGT')
        consT = request.form.get('consT')
        posion_start = request.form.get('posion_start')
        posion_end = request.form.get('posion_end')
        GCmin = request.form.get('GCmin')
        GCmax = request.form.get('GCmax')
        first = request.form.get('first')
        tC = request.form.get('tC')
        tA = request.form.get('tA')
        tG = request.form.get('tG')
        # database = request.form.get('database')
        file = request.files['file']

        blastCheck = request.form.get('blastCheck')
        sevenbp = request.form.get('sevenbp')
        # 在这里进行其他处理...

        # 删除所有的空格和数字的sequence
        sequence = sequence.replace(' ', '')
        sequence = ''.join([i for i in sequence if not i.isdigit()])
        seq = sequence.upper()  # 转为大写
        seq = seq.replace('\n', '').replace('\r', '')  # 移除换行符
        sequence_length = len(seq)  # 序列长度
        db_info = 0

        # 预测靶基因的二级结构8
        with open('input.txt', 'w') as f:
            f.write(seq)
        p = subprocess.Popen('.\\tools\\ViennaRNA\\RNAfold.exe --noPS < input.txt > secondary.txt',
                             shell=True)
        p.wait()

        # 处理文件上传，并建库

        if file:
            # 在这里，你可以处理文件，例如保存到服务器
            # file.save('./uploads/file.fq')
            filename = 'file.fq'
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            # 建库
            db_path = Path(f'./uploads/file.fq')
            makeblastdb_path = Path('./tools/ncbi_blast/bin/makeblastdb.exe')
            out_path = Path('uploads/blastdb/file/file')
            command = f'{makeblastdb_path} -in {db_path} -dbtype nucl -parse_seqids -out {out_path} -title "file"'
            os.system(command)
            db_info = 1  # 用于告知建库，result.html

        # 检查序列是否正确
        def check_sequence(sequence):
            valid_chars = set('CGTAU')
            for char in sequence:
                if char not in valid_chars:
                    return False
            return True


        # 排除连续碱基的序列
        def check_seq_gswn(seq):
            # 排除包含回文序列的序列
            # 排除已知诱导免疫反应的序列GTCCTTCAA、TGTGT和CTGAATT
            seq = seq.replace('U', 'T')

            if consGC == '' and consAT == '' and consT == '' and consGT == '':
                return not (seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)
            elif consGC == '' and consT == '' and consGT == '':
                return not (re.search(f'[AT]{{{consAT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)
            elif consGC == '' and consAT == '' and consGT == '':
                return not (re.search(f'[T]{{{consT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)
            elif consGC == '' and consAT == '' and consT == '':
                return not (re.search(f'[GT]{{{consGT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)



            elif consGC == '' and consAT == '':
                return not (re.search(f'[GT]{{{consGT},}}', seq) or
                            re.search(f'[T]{{{consT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)
            elif consGC == '' and consT == '':
                return not (re.search(f'[AT]{{{consAT},}}', seq) or
                            re.search(f'[GT]{{{consGT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)
            elif consGC == '' and consGT == '':
                return not (re.search(f'[AT]{{{consAT},}}', seq) or
                            re.search(f'[GT]{{{consGT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)

            elif consGC == '':
                return not (re.search(f'[GT]{{{consGT},}}', seq) or
                            re.search(f'[AT]{{{consAT},}}', seq) or
                            re.search(f'[T]{{{consT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)

            elif consAT == '' and consT == '' and consGT == '':
                return not (re.search(f'[GC]{{{consGC},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)
            elif consAT == '' and consT == '':
                return not (re.search(f'[GC]{{{consGC},}}', seq) or
                            re.search(f'[GT]{{{consGT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)
            elif consAT == '' and consGT == '':
                return not (re.search(f'[GC]{{{consGC},}}', seq) or
                            re.search(f'[T]{{{consT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)
            elif consAT == '':
                return not (re.search(f'[GT]{{{consGT},}}', seq) or
                            re.search(f'[GC]{{{consGC},}}', seq) or
                            re.search(f'[T]{{{consT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)

            elif consT == '' and consGT == '':
                return not (re.search(f'[GC]{{{consGC},}}', seq) or
                            re.search(f'[AT]{{{consAT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)
            elif consT == '':
                return not (re.search(f'[GT]{{{consGT},}}', seq) or
                            re.search(f'[GC]{{{consGC},}}', seq) or
                            re.search(f'[AT]{{{consAT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)

            elif consGT == '':
                return not (re.search(f'[GC]{{{consGC},}}', seq) or
                            re.search(f'[AT]{{{consAT},}}', seq) or
                            re.search(f'[T]{{{consT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)

            else:
                return not (re.search(f'[GT]{{{consGT},}}', seq) or
                            re.search(f'[AT]{{{consGC},}}', seq) or
                            re.search(f'[T]{{{consAT},}}', seq) or
                            re.search(f'[GC]{{{consT},}}', seq) or
                            seq == seq[::-1] or
                            'GTCCTTCAA' in seq or
                            'TGTGT' in seq or
                            'CTGAATT' in seq)
            # # 排除连续4个或更多的T的序列
            # # 排除连续6个或更多的G和C或G和T的序列
            # return not (re.search('[TU]{4,}', seq) or
            #             re.search('[GC]{6,}', seq) or
            #             re.search('[GTU]{6,}', seq) or
            #             seq == seq[::-1] or
            #             'GTCCTTCAA' in seq or
            #             'TGTGT' in seq or
            #             'CTGAATT' in seq)
        # 字符长度保持一样长，填充短的字符串
        def pad_strings(str1, str2):
            len1 = len(str1)
            len2 = len(str2)
            if len1 > len2:
                str2 += ' ' * (len1 - len2)
            elif len2 > len1:
                str1 += ' ' * (len2 - len1)
            return str1, str2

        # 互补碱基
        def complement(seq):
            complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
            return "".join(complement[base] for base in seq)

        # 打分函数,总共
        def dafen(seqTU):
            length = abs(int(taille_motif))
            senseq = seqTU[-length:]
            antiseq = complement(seqTU)  # antisense strand 3'-5'
            antiseq = antiseq[:length]
            antiseq = antiseq[::-1]  # antisense strand 5'-3'

            score_list = []  # 得分详情，最后一位是总得分

            # （1）反义链第二至第七个核苷酸中的GC含量在33%到67%( 这段长度6的序列里GC数量为2-4)，
            # 反义链第八至第十八个核苷酸中的GC含量应为33%到55%(GC数量为4-6)
            antiGC27 = antiseq[1: 7].count('G') + antiseq[1: 7].count('C')
            antiGC818 = antiseq[7: 18].count('G') + antiseq[7: 18].count('C')
            if 2 < antiGC27 < 4 and 4 < antiGC818 < 6:
                score_list.append(1)
            else:
                score_list.append(0)

            # (2)antisense GC重复次数不超过1次，AT重复次数不超过2次
            gc_repeat = antiseq.count('GC') - 1
            at_repeat = antiseq.count('AT') - 1
            if gc_repeat < 2 and at_repeat < 3:
                score_list.append(1)
            else:
                score_list.append(0)

            # (3)碱基A连续重复最大次数不超过3次，
            # 碱基T连续重复最大次数不超过3次，或者U
            # 碱基G连续重复最大次数不超过2次，
            # 碱基C连续重复最大次数不超过2次
            match_A = re.search('[A]{4,}', antiseq)
            match_U = re.search('[U]{4,}', antiseq)
            match_G = re.search('[G]{3,}', antiseq)
            match_C = re.search('[C]{3,}', antiseq)

            if match_A is None and \
                    match_U is None and \
                    match_G is None and \
                    match_C is None:

                score_list.append(1)
            else:
                score_list.append(0)
            # （4）TT overhangs in 3’-end
            # 计算悬挂
            sense_overhangs = senseq[-2:]
            antiseq__overhangs = antiseq_reversed[-2:]
            if sense_overhangs == 'TT' and antiseq__overhangs == 'TT':

                score_list.append(1)
            else:
                score_list.append(0)


            # 范围 18-29
            if 17 < length < 30:
                # (5)反义链3′端碱基 - --26nt的3’端碱基为A或T；----1其他长度就直接 + 1
                if length == 26:
                    if antiseq[-1] in ['A', 'U']:

                        score_list.append(1)
                    else:
                        score_list.append(0)
                else:

                    score_list.append(1)

                # (6)选取29长度时，加上第19位置上碱基为A或T的特征。----1        其他长度就直接+1
                if length == 29:
                    if antiseq[18] in ['A', 'U']:

                        score_list.append(1)
                    else:
                        score_list.append(0)
                else:

                    score_list.append(1)

            # 范围 30-40
            elif 29 < length < 41:
                # （5）反义链3′端碱基---33nt的3’端碱基为T。
                if length == 33:
                    if antiseq[-1] in ['U']:

                        score_list.append(1)
                    else:
                        score_list.append(0)
                else:

                    score_list.append(1)


                # (6)第21个碱基位为A，---1
                if antiseq[20] == 'A':

                    score_list.append(1)
                else:

                    score_list.append(0)

                # (7)第25个碱基位为T或C，---1
                if antiseq[24] in ['U', 'C']:

                    score_list.append(1)
                else:

                    score_list.append(0)

                # (8)第30个碱基位为T或G，---1
                if antiseq[29] in ['U', 'G']:

                    score_list.append(1)
                else:

                    score_list.append(0)

                if length > 32:
                    # (9)第33个碱基位为T。---1
                    if antiseq[32] in ['U']:

                        score_list.append(1)
                    else:

                        score_list.append(0)
                else:
                    score_list.append(0)
                # (10)对于设计31、32长度时，再附加上第10碱基位为G的特征--1
                if length == 32:
                    if antiseq[-10] in ['G']:

                        score_list.append(1)
                    else:

                        score_list.append(0)
                else:
                    score_list.append(1)

            return score_list

        # （7）antisense的5‘端->3'端的前八个位置，对应的靶基因的位置的位置，要是都是点4分，要是前八个是点就2分
        def targetSec(startS):
            startS = int(startS)  # 起始位置0
            score = 0
            endS = int(taille_motif) - 1 + startS  # 末尾位置
            with open('secondary.txt', 'r') as f:
                lines = f.readlines()[1]
                second_line = lines[startS - 1:endS]  # 整行，不包括能量，antisense链对应的靶基因区域的二级结构
                # all_count = second_line.count('.')  # 全部点数有多少
                eight_second = second_line[-8:]  # antisensed的5‘-3’方向的对应的靶基因二级结构的前八个位置
                count_dian = eight_second.count('.')  # 前八个点数量
                if count_dian == 8:
                    score += 2
            return score

        # 靶位点的二级结构
        def targetStruc(startS):
            startS = int(startS)  # 起始位置
            score = 0
            endS = int(taille_motif) - 1 + startS + 2 # 末尾位置
            with open('secondary.txt', 'r') as f:
                lines = f.readlines()[1]
                second_line = lines[startS - 1:endS]  # 整行，不包括能量，antisense链对应的靶基因区域的二级结构
            return second_line





        # Check if the application is running as a PyInstaller bundle
        if getattr(sys, 'frozen', False):
            application_path = sys._MEIPASS
        else:
            application_path = Path.cwd()

        blastn_path = Path(application_path, 'tools/ncbi_blast/bin/blastn.exe')
        # 设置 blastdb 目录的相对路径
        db_path = Path(application_path, 'uploads/blastdb/file/file')

        # 数据库blastn
        def run_blastn(blast_seq_file):
            word_size = 13
            evalue = 1
            command = f'{blastn_path} -query {blast_seq_file} -out blast_output.txt -db {db_path} -word_size {word_size} -evalue {evalue} -task blastn-short'
            os.system(command)

        # 统计最长连续相同的碱基对长度
        def longest_match(query, sbjct):
            matches = [q == s for q, s in zip(query, sbjct)]
            count = max_count = 0
            for match in matches:
                if match:
                    count += 1
                    max_count = max(max_count, count)
                else:
                    count = 0
            return max_count



        if check_sequence(seq) == False:
            message_error = "Error! The sequence contains unauthorized characters."
            return render_template('result.html', message_error=message_error)
        # 计算序列的长度，并进行判断
        elif sequence_length == 0:  # 添加错误检查，确保sequence变量不为None或空字符串
            message_error = "Enter a sequence, please"
            return render_template('result.html', message_error=message_error)
        elif sequence_length < 22:
            message_error = "The sequence is too short"
            return render_template('result.html', message_error=message_error)
        else:
            # 添加错误检查，确保posion_start和posion_end都存在并能够正确转换为整数
            if posion_start == "" and posion_end == "":
                sequence_length = len(seq)
            elif posion_start == "" and posion_end != "":
                try:
                    posion_end = abs(int(posion_end))
                    if posion_end > len(seq):
                        message_error = "Error, position values are not valid"
                        return render_template('result.html', message_error=message_error)
                    elif posion_end < 22:
                        message_error = "Error, position values are not valid"
                        return render_template('result.html', message_error=message_error)
                    else:
                        sequence_length = posion_end
                        seq = seq[:posion_end]
                except ValueError:
                    message_error = "Error, position values are not valid"
                    return render_template('result.html', message_error=message_error)
            elif posion_start != "" and posion_end == "":
                try:
                    posion_start = abs(int(posion_start))
                    if posion_start > len(seq):
                        message_error = "Error, position values are not valid"
                        return render_template('result.html', message_error=message_error)
                    elif (len(seq) - posion_start + 1) < 22:
                        message_error = "Error, position values are not valide"
                        return render_template('result.html', message_error=message_error)
                    else:
                        sequence_length = len(seq) - posion_start + 1
                        seq = seq[posion_start:]
                except ValueError:
                    message_error = "Error, position values are not valide"
                    return render_template('result.html', message_error=message_error)

            else:
                try:
                    posion_start = abs(int(posion_start))
                    posion_end = abs(int(posion_end))
                    if posion_start > len(seq):
                        message_error = "Error, position values are not valid"
                        return render_template('result.html', message_error=message_error)
                    elif (len(seq) - posion_start) < 22:
                        message_error = "Error, position values are not valide"
                        return render_template('result.html', message_error=message_error)
                    elif posion_end > len(seq):
                        message_error = "Error, position values are not valid"
                        return render_template('result.html', message_error=message_error)
                    elif posion_end - posion_start < 22:
                        message_error = "Error, position values are not valid"
                        return render_template('result.html', message_error=message_error)
                    else:
                        sequence_length = posion_end - posion_start + 1
                        seq = seq[posion_start:posion_end]
                except ValueError:
                    message_error = "Error, position values are not valide"
                    return render_template('result.html', message_error=message_error)

        # GC含量的判断
        if GCmin == "":
            message_error = "Error, GC% values are not valid"
            return render_template('result.html', message_error=message_error)
        elif GCmax == "":
            message_error = "Error, GC% values are not valid"
            return render_template('result.html', message_error=message_error)
        else:
            try:
                GCmin = abs(int(GCmin))
                GCmax = abs(int(GCmax))
                if GCmin > 100:
                    message_error = "Error, GC% values are not valid"
                    return render_template('result.html', message_error=message_error)
                elif GCmax > 100:
                    message_error = "Error, GC% values are not valid"
                    return render_template('result.html', message_error=message_error)
                elif GCmin >= GCmax:
                    message_error = "Error, GC% values are not valid"
                    return render_template('result.html', message_error=message_error)
            except ValueError:
                message_error = "Error, GC% values are not valid"
                return render_template('result.html', message_error=message_error)
        # motif size进行判断 18-40,"Error! Please, enter a valid motif size" posion_start,posion_end

        if taille_motif:
            try:
                taille_motif = abs(int(taille_motif))
                if taille_motif < 18:
                    message_error = "Error! Please, enter a valid motif size"
                    return render_template('result.html', message_error=message_error)
                elif taille_motif > 40:
                    message_error = "Error! Please, enter a valid motif size"
                    return render_template('result.html', message_error=message_error)
            except ValueError:
                message_error = "Error! Please, enter a valid motif size"
                return render_template('result.html', message_error=message_error)
        else:
            message_error = "Error! Please, enter a valid motif size"
            return render_template('result.html', message_error=message_error)

        ''' 开始进行计算 '''
        # 图形滑动窗口法
        length = abs(int(taille_motif))
        winlen = length + 2  # 图形滑动窗口长度，需要+2，因为有悬挂

        """
        长度18-29       (8分）
        (1)反义链第二至第七个核苷酸中的GC含量在33%到67%( 这段长度6的序列里GC数量为2-4)，反义链第八至第十八个核苷酸中的GC含量应为33%到55%(GC数量为4-6)----2
        (2)antisense GC重复次数不超过1次，AT重复次数不超过2次-----1
        (3)碱基A连续重复最大次数不超过3次，碱基T连续重复最大次数不超过3次，或者U碱基G连续重复最大次数不超过2次，碱基C连续重复最大次数不超过2次-------2
        （4）TT overhangs in 3’-end-----1
        (5)反义链3′端碱基---26nt的3’端碱基为A或T；----1        其他长度就直接+1
        (6)选取29长度时，加上第19位置上碱基为A或T的特征。----1        其他长度就直接+1
        
        
        长度30-40      (12分）
        (1)反义链第二至第七个核苷酸中的GC含量在33%到67%( 这段长度6的序列里GC数量为2-4)，反义链第八至第十八个核苷酸中的GC含量应为33%到55%(GC数量为4-6)----2
        (2)antisense GC重复次数不超过1次，AT重复次数不超过2次-----1
        (3)碱基A连续重复最大次数不超过3次，碱基T连续重复最大次数不超过3次，或者U碱基G连续重复最大次数不超过2次，碱基C连续重复最大次数不超过2次-------2
        (4)TT overhangs in 3’-end-----1
        (5)反义链3′端碱基---33nt的3’端碱基为T。----1
        (6)第21个碱基位为A，---1
        (7)第25个碱基位为T或C，---1
        (8)第30个碱基位为T或G，---1
        (9)第33个碱基位为T。---1
        (10)对于设计31、32长度时，再附加上第10碱基位为G的特征--1        其他长度就直接+1
        """
        score_columns = []  # 得分情况，列名
        if 17 < length < 30:
            for i in range(6+1):  # 加上一个二级结构的判断
                score_columns.append(str(i+1))
        elif 29 < length < 41:
            for i in range(10+1):  # 加上一个二级结构的判断
                score_columns.append(str(i+1))

        # # 列名
        # table_columns = ['target sequence <br>' + str(taille_motif) + 'nt target + 2nt overhang',
        #                  'target sequence <br> secondary structure', 'target position',
        #                  str(taille_motif) + 'nt passenger(5′→3′)',
        #                  str(taille_motif) + 'nt guide(5′→3′)', 'GC%', 'Score', ]  # 表格的列名
        # 列名
        table_columns = ['target sequence <br>' + str(taille_motif) + 'nt target + 2nt overhang',
                         'target sequence <br> secondary structure', 'target position',
                         str(taille_motif) + 'nt passenger(5′→3′)',
                         str(taille_motif) + 'nt guide(5′→3′)', 'GC%', 'Score', ]  # 表格的列名
        table_columns.extend(score_columns)
        # 使用Markup函数将每个字符串转换为Markup对象
        table_columns = [Markup(column) for column in table_columns]
        seq_table = []  # 存放符合条件的序列

        # 注意Restrict the search to the region from position，需要设置起始点和结束位置

        for a in range(sequence_length):
            index_left = a
            index_right = index_left + winlen
            if index_right > sequence_length - 1:  # 如果窗口的右边界超过了序列的长度，就结束循环
                break
            elif index_right < (winlen - 1):  # 如果窗口的右边界小于窗口长度减一，也结束循环
                break
            seq1 = seq[index_left:index_right]
            seqTU = seq1.replace('T', 'U')  # 把T换成U

            senseq = seqTU[-length:]

            antiseq = complement(seqTU)  # antisense strand 3'-5'
            antiseq = antiseq[:length]
            antiseq_reversed = antiseq[::-1]  # antisense strand 5'-3'

            if check_seq_gswn(seq1):
                countGC = (antiseq.count('G') + antiseq.count('C')) / winlen * 100
                formatted_GC = "{:.2f}".format(countGC)

                # 符合条件的留下----antisense的GC含量在GCmin和GCmax之间
                if int(GCmin) <= float(formatted_GC) <= int(GCmax):  # GC含量

                    # 符合条件的靶基因的序列的第一个核苷酸是first
                    if first == '':
                        score_list = dafen(seqTU)
                        score = sum(score_list)
                        # 使用 map() 函数和 str() 函数将每个元素转换为字符串
                        score_list = list(map(str, score_list))

                        if posion_start == '':
                            secondStr = targetStruc(index_left + 1)
                            secondScore = targetSec(index_left + 1)
                            score = score + secondScore
                            if 17 < length < 30:
                                if score > 4:
                                    seq_table.append(
                                        dict(zip(table_columns,
                                                 [seq1, secondStr, str(index_left + 1) + '-' + str(index_right),
                                                  senseq, antiseq_reversed, formatted_GC, score]
                                                 + score_list + [str(secondScore)])))

                            # 范围 30-40
                            elif 29 < length < 41:
                                if score > 8:
                                    seq_table.append(
                                        dict(zip(table_columns,
                                                 [seq1, secondStr, str(index_left + 1) + '-' + str(index_right),
                                                  senseq, antiseq_reversed, formatted_GC, score]
                                                 + score_list + [str(secondScore)])))
                        else:
                            secondStr = targetStruc(posion_start)
                            secondScore = targetSec(posion_start)
                            score = score + secondScore
                            if 17 < length < 30:
                                if score > 5:
                                    seq_table.append(
                                        dict(zip(table_columns,
                                                 [seq1, secondStr, str(posion_start+1) + '-' + str(posion_start + winlen),
                                                  senseq, antiseq_reversed, formatted_GC, score]
                                                 + score_list + [str(secondScore)])))


                            # 范围 30-40
                            elif 29 < length < 41:
                                if score > 8:
                                    seq_table.append(
                                        dict(zip(table_columns,
                                                 [seq1, secondStr, str(posion_start+1) + '-' + str(posion_start + winlen),
                                                  senseq, antiseq_reversed, formatted_GC, score]
                                                 + score_list + [str(secondScore)])))


                    elif first == 'A or G':
                        if seq1[0] == 'A' and seq1[-1] not in ['G', 'C']:
                            score_list = dafen(seqTU)
                            score = sum(score_list)
                            # 使用 map() 函数和 str() 函数将每个元素转换为字符串
                            score_list = list(map(str, score_list))

                            if posion_start == '':
                                secondStr = targetStruc(index_left + 1)
                                secondScore = targetSec(index_left + 1)
                                score = score + secondScore
                                if 17 < length < 30:
                                    if score > 4:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr, str(index_left + 1) + '-' + str(index_right),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))

                                # 范围 30-40
                                elif 29 < length < 41:
                                    if score > 8:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr, str(index_left + 1) + '-' + str(index_right),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))
                            else:
                                secondStr = targetStruc(posion_start)
                                secondScore = targetSec(posion_start)
                                score = score + secondScore
                                if 17 < length < 30:
                                    if score > 5:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr,
                                                      str(posion_start + 1) + '-' + str(posion_start + winlen),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))


                                # 范围 30-40
                                elif 29 < length < 41:
                                    if score > 8:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr,
                                                      str(posion_start + 1) + '-' + str(posion_start + winlen),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))
                        elif seq1[0] == 'G':
                            score_list = dafen(seqTU)
                            score = sum(score_list)
                            # 使用 map() 函数和 str() 函数将每个元素转换为字符串
                            score_list = list(map(str, score_list))
                            # if db_info == 1:
                            #     score = dafen(seqTU)
                            # else:
                            #     score = dafen(seqTU) + 4
                            if posion_start == '':
                                secondStr = targetStruc(index_left + 1)
                                secondScore = targetSec(index_left + 1)
                                score = score + secondScore
                                if 17 < length < 30:
                                    if score > 4:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr, str(index_left + 1) + '-' + str(index_right),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))

                                # 范围 30-40
                                elif 29 < length < 41:
                                    if score > 8:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr, str(index_left + 1) + '-' + str(index_right),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))
                            else:
                                secondStr = targetStruc(posion_start)
                                secondScore = targetSec(posion_start)
                                score = score + secondScore
                                if 17 < length < 30:
                                    if score > 5:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr,
                                                      str(posion_start + 1) + '-' + str(posion_start + winlen),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))


                                # 范围 30-40
                                elif 29 < length < 41:
                                    if score > 8:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr,
                                                      str(posion_start + 1) + '-' + str(posion_start + winlen),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))
                    elif first == 'A':
                        if seq1[0] == first and seq1[-1] not in ['G', 'C']:
                            score_list = dafen(seqTU)
                            score = sum(score_list)
                            # 使用 map() 函数和 str() 函数将每个元素转换为字符串
                            score_list = list(map(str, score_list))
                            # if db_info == 1:
                            #     score = dafen(seqTU)
                            # else:
                            #     score = dafen(seqTU) + 4
                            if posion_start == '':
                                secondStr = targetStruc(index_left + 1)
                                secondScore = targetSec(index_left + 1)
                                score = score + secondScore
                                if 17 < length < 30:
                                    if score > 4:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr, str(index_left + 1) + '-' + str(index_right),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))

                                # 范围 30-40
                                elif 29 < length < 41:
                                    if score > 8:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr, str(index_left + 1) + '-' + str(index_right),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))
                            else:
                                secondStr = targetStruc(posion_start)
                                secondScore = targetSec(posion_start)
                                score = score + secondScore
                                if 17 < length < 30:
                                    if score > 5:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr,
                                                      str(posion_start + 1) + '-' + str(posion_start + winlen),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))


                                # 范围 30-40
                                elif 29 < length < 41:
                                    if score > 8:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr,
                                                      str(posion_start + 1) + '-' + str(posion_start + winlen),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))
                    elif first == 'G':
                        if seq1[0] == first:
                            score_list = dafen(seqTU)
                            score = sum(score_list)
                            # 使用 map() 函数和 str() 函数将每个元素转换为字符串
                            score_list = list(map(str, score_list))
                            # if db_info == 1:
                            #     score = dafen(seqTU)
                            # else:
                            #     score = dafen(seqTU) + 4
                            if posion_start == '':
                                secondStr = targetStruc(index_left + 1)
                                secondScore = targetSec(index_left + 1)
                                score = score + secondScore
                                if 17 < length < 30:
                                    if score > 4:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr, str(index_left + 1) + '-' + str(index_right),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))

                                # 范围 30-40
                                elif 29 < length < 41:
                                    if score > 8:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr, str(index_left + 1) + '-' + str(index_right),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))
                            else:
                                secondStr = targetStruc(posion_start)
                                secondScore = targetSec(posion_start)
                                score = score + secondScore
                                if 17 < length < 30:
                                    if score > 5:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr,
                                                      str(posion_start + 1) + '-' + str(posion_start + winlen),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))


                                # 范围 30-40
                                elif 29 < length < 41:
                                    if score > 8:
                                        seq_table.append(
                                            dict(zip(table_columns,
                                                     [seq1, secondStr,
                                                      str(posion_start + 1) + '-' + str(posion_start + winlen),
                                                      senseq, antiseq_reversed, formatted_GC, score]
                                                     + score_list + [str(secondScore)])))
            if posion_start != '':
                posion_start += 1




        candidate = len(seq_table)
        if candidate == 0:
            message_error = "Empty! No matching siRNA candidate"
            return render_template('result.html', message_error=message_error)


        # blast 分析
        seq_table_blast = seq_table
        seq_off_table = []  # 存放未脱靶的序列

        # seq_off_blast = []  # 存放未脱靶序列的blast结果

        # 数据库建立好，进行blast分析
        if db_info == 1:
            table_columns.append('off-target')
            table_columns.append('show')
            for s in seq_table_blast:
                seq = s['target sequence <br>' + str(taille_motif) + 'nt target + 2nt overhang']
                seqTU = seq.replace('T', 'U')  # 把T换成U
                antiseq = complement(seqTU)  # antisense strand 3'-5'
                antiseq = antiseq[:length]
                antiseq = antiseq[::-1]  # antisense strand 5'-3'


                with open('blast_seq.fastq', 'w') as f:
                    f.write(antiseq)
                run_blastn('blast_seq.fastq')

                guidecount = 0
                with open('blast_output.txt', 'r') as f:
                    content = f.read()  # 读取文件内容

                    bpmatchNum = []  # 存放连续核苷酸相同的碱基数量
                    queryNum = []  # 存放blast-antisense覆盖率
                    antimatchNum = []  # 存放antisense_blast的匹配数目

                    query_list = []  # 存放查询的序列
                    su_list = []  # 存放匹配程度的线条
                    sbjct_list = []  # 存放比对到的序列
                    seqid = ''  # “>” 后面的部分是序列的 ID

                    # 初始化同源序列的计数
                    blastCount = 0
                    query = sbjct = ''

                    # 检查是否存在'Strand=Plus/Minus'
                    if 'Strand=Plus/Minus' in content:

                        for line in content.split('\n'):

                            # 如果行以 ">" 开头，表示一个新的同源序列
                            if line.startswith(">"):

                                seqid = line[1:].split()[0]# 提取序列ID，忽略 ">" 并取第一个空格分隔的部分

                                blastCount += 1

                            # 检查该行是否包含所需的内容
                            elif line.startswith('Query  '):
                                print(line)
                                # query = line.split()[2]#只获取序列
                                # query_list.append(query)
                                guidecount += 1
                                guide = 'guide' + str(guidecount)
                                seqid, guide = pad_strings(seqid, guide)

                                line = line.replace('Query', guide)
                                query_list.append(line)

                            elif "||" in line:  # 使用正则表达式匹配匹配程度的线条
                                print(line)
                                line = ' ' * (len(seqid) - 5) + line

                                su_list.append(line)

                            elif line.startswith('Sbjct  '):
                                print(line)
                                # sbjct = line.split()[2] #只获取序列
                                # sbjct_list.append(sbjct)

                                # 替换为ID

                                line = line.replace('Sbjct', seqid)
                                sbjct_list.append(line)

                                bpmatch = longest_match(query, sbjct)
                                bpmatchNum.append(bpmatch)

                        if bpmatchNum:
                            # 候选siRNAs中有超过13个连续核苷酸与另一个mRNA的序列相同的被丢弃，同源序列超过5条
                            # if all(num < 14 for num in bpmatchNum) and blastCount < 6:
                            # 候选siRNAs同源序列超过5条就丢弃
                            if blastCount < 6:
                                s['off-target'] = blastCount

                                # 使用zip函数将三个列表中相同索引的元素组合在一起
                                seq_off_blast = ['<br>'.join(item) for item in zip(query_list, su_list, sbjct_list)]
                                # s['show'] = seq_off_blast[0]
                                s['show'] = '<br>'.join(seq_off_blast)

                                seq_off_table.append(s)



                    elif 'No hits found' in content:
                        print('No hits found')

                        query_list.append(' ')
                        su_list.append(' ')
                        sbjct_list.append(' ')

                        s['off-target'] = blastCount
                        # 使用zip函数将三个列表中相同索引的元素组合在一起
                        seq_off_blast = ['<br>'.join(item) for item in zip(query_list, su_list, sbjct_list)]
                        s['show'] = '<br>'.join(seq_off_blast)

                        seq_off_table.append(s)



            os.remove('blast_seq.fastq')
            os.remove('blast_output.txt')

        os.remove('input.txt')
        os.remove('secondary.txt')

        # # 使用Markup函数将每个字符串转换为Markup对象
        # seq_off_blast = [Markup(content) for content in seq_off_blast]

        off_candidate = len(seq_off_table)
        if db_info == 1 and off_candidate == 0:
            message_error = "SiRNA candidate is all off target"
            return render_template('result.html', message_error=message_error)
        seq_table = sorted(seq_table, key=lambda x: x['Score'], reverse=True)
        seq_off_table = sorted(seq_off_table, key=lambda x: x['Score'], reverse=True)
        """将结果传递给result模板"""

        # 将序列长度，siRNA靶序列GC含量传递给模板
        return render_template('result.html', nom=nom, sequence_length=sequence_length,
                               GCmin=GCmin, GCmax=GCmax, first=first, taille_motif=taille_motif,
                               message_tC=tC, message_tA=tA, message_tG=tG, score_columns=score_columns,
                               table_columns=table_columns, seq_table=seq_table, candidate=candidate,
                               seq_off_table=seq_off_table, off_candidate=off_candidate, db_info=db_info)


    else:
        # 这是一个GET请求，渲染表单页面
        return render_template('siRNA_designer.html')


if __name__ == "__main__":
    app.run()
