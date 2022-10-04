# by Francisco Ortega, June 2008.
# Seismological Laboratory
# California Institute of Technology
# add more comments here

class LogFile:
   """
   Defines an object that manages log files, to be used in any python script.

   usage:
   from FileClasses import LogFile
   log = LogFile('logFileName')
   log.clearFile() # dont use if you want to append the logs to an existing file.
   log.addLine(' message to log') # can be a multiline string.

   """
   def __init__(self, filename):
      self.filename = filename

   def getFileName(self):
      return self.filename

   def setFileName(self,filename):
      self.filename = filename

   def clearFile(self):
      import time
      with open(self.filename,'w') as fout:
          fout.write(time.ctime() + ":\n-> Log file initialization time\n")
      print(time.ctime() + ":\n-> Log file initialization time\n")
      

   def addLine(self,textline):
      import time
      with open(self.filename,'a') as fout:
          fout.write(time.ctime() + ":\n")
          fout.write("-> " + textline+"\n")
      print(time.ctime() + ":")
      print("-> " + textline)

