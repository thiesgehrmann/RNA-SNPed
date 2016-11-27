import os;
import sys;
import shlex, subprocess;

def run_cmd(cmd, bg = False, stdin = None, stdout = None, stderr = None):
    print '[RUN CMD] Executing: %s' % cmd
    sys.stdout.flush()
    p = subprocess.Popen(shlex.split(cmd), stdin=stdin, stdout=stdout, stderr=stderr);
    if bg :
        return p
    else:
        (pid, r) = os.waitpid(p.pid, 0);
        return r;

def run_cmd_fail(cmd, stdin = None, stdout = None, stderr = None):
    print '[RUN CMD] Executing: %s' % cmd
    sys.stdout.flush()
    p = subprocess.Popen(shlex.split(cmd), stdin=stdin, stdout=stdout, stderr=stderr);
    (pid, r) = os.waitpid(p.pid, 0);
    if r != 0:
      raise RuntimeError, "Command failed: " + str(res) + "\n" + cmd;
      sys.exit(r);
    return r;

def run_shell(cmd, bg = False):
    print '[RUN SHELL] Executing: %s' % cmd
    sys.stdout.flush()
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    if bg :
        return p
    else:
        (pid, r) = os.waitpid(p.pid, 0);
        return r;

def getCommandOutput(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    if p.returncode != 0:
        raise RuntimeError, "Command failed: " + str(res) + "\n" + cmd;
    return res
