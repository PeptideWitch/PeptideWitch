from flask import Flask, render_template, request, send_file
import os
import datetime
import shutil

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def home(name=None):
    if request.method == "GET":
        if os.path.isfile("Output.zip"):
            print("File Found")
            os.remove("Output.zip")
    return render_template('home.html', name=name)

@app.route('/upload', methods=['GET', 'POST'])
def upload():
    dayandTime = str(datetime.datetime.now().strftime('%d-%m-%Y_%H-%M-%S'))
    direc = os.path.join(os.getcwd(), 'uploads' , dayandTime)
    try:
        os.mkdir(direc)
    except FileExistsError:
        print("Directory exists")
    print(direc)
    if request.method == "POST":
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
        os.system("python PepWitch1.3.Online.py " + str(minspc) + " " + str(disregard) + " "  + str(pval) + " " + str(replicates) + " " + str(spcfrac) + " " + str(direc) + " "+ str(names))
        shutil.rmtree(direc)
        return send_file("Output.zip", as_attachment=True)

if __name__ == "__main__":
    app.run()