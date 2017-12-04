from flask import Flask, render_template, request, send_file, url_for
import os
import datetime
import shutil
import sys

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def home(name=None):
    if request.method == "GET":
        if os.path.isfile("Output.zip"):
            os.remove("Output.zip")
            print("File Erased")
        return render_template('home2.html', name=name)
    if request.method == "POST":
        print("Here" + sys.version)
        DayandTime = str(datetime.datetime.now().strftime('%d-%m-%Y_%H-%M-%S'))
        direc = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'uploads' , DayandTime)
        try:
            os.mkdir(direc)
        except FileExistsError:
            print("Directory exists")
        print(direc)
        uploaded_files = request.files.getlist("file[]")
        names = request.form["names"]
        minspc = request.form["minspc"]
        disregard = request.form["disregard"]
        pval = request.form["pval"]
        replicates = request.form["replicate"]
        spcfrac = request.form["spcfrac"]
        for file in uploaded_files:
            newpath = os.path.join(direc, file.filename)
            file.save(newpath)
            file.close()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        os.system("python -V")
        os.system("/home/PeptideWitch/.virtualenvs/my-virtualenv/bin/python3.6 PepWitch1.3.Online.py " + str(minspc) + " " + str(disregard) + " "  + str(pval) + " " + str(replicates) + " " + str(spcfrac) + " " + str(direc) + " "+ str(names))
        shutil.rmtree(direc)
        return send_file(os.path.dirname(os.path.realpath(__file__)) + "/Output.zip", as_attachment=True)

if __name__ == "__main__":
    app.run()
