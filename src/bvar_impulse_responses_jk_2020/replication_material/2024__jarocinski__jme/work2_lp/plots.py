# Create a latex file that includes all pdf figures from folder 'pdfdir'
import os
import glob

tex_code = r"""\documentclass[11pt]{article}
\usepackage{amssymb}
\usepackage{amsfonts}
\usepackage{amsmath}
\usepackage{graphicx}
\usepackage[round]{natbib}
\usepackage{verbatim}
\usepackage{fancyhdr}
\usepackage{booktabs}
\usepackage[hmargin=0.5cm,vmargin=2cm]{geometry}
\begin{document}
\pagestyle{fancy}
\fancyhead[C]{}
\fancyfoot[C]{\thepage}
\noindent"""

pdfdir = 'out'

outtex_fname = pdfdir + ".tex"
if os.path.isfile(outtex_fname):
	os.remove(outtex_fname)
    
ncol = 4
width = 0.22
cc = 0
for pdfpath in glob.glob(pdfdir + "/*.pdf"):
    fname = pdfpath.replace("\\","/")
    tex_code += '\\includegraphics[width=' + str(width) + '\\textwidth]{' + fname + '}\n'
    cc += 1
    if cc == ncol:
        tex_code += '\\\\[3ex]\n'
        cc = 0

tex_code += '\end{document}'

tex_file = open(outtex_fname,'w')
tex_file.write(tex_code)
tex_file.close()

os.system('texify -c -p ' + outtex_fname)

