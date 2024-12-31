import paramiko
import pandas as pd
import xml.etree.ElementTree as ET
from getpass import getpass
import time

class SSHConnection:
    def __init__(self, hostname, username, key_path=None, password=None):
        self.hostname = hostname
        self.username = username
        self.key_path = key_path
        self.password = password
        self.client = None

    def __enter__(self):
        self.client = paramiko.SSHClient()
        self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        try:
            # First attempt: Try SSH key authentication
            if not self.key_path:
                self.key_path = f"/home/{self.username}/.ssh/id_rsa"
            
            try:
                self.client.connect(self.hostname, 
                                  username=self.username,
                                  key_filename=self.key_path)
            except (paramiko.ssh_exception.AuthenticationException, FileNotFoundError):
                # Second attempt: Fall back to password authentication
                if not self.password:
                    self.password = getpass(f"Enter password for {self.username}@{self.hostname}: ")
                self.client.connect(self.hostname, 
                                  username=self.username,
                                  password=self.password)
            
            return self
        except Exception as e:
            print(f"Failed to connect to {self.hostname}: {str(e)}")
            return None
    
    def execute_commands(self, commands):
        """Execute multiple commands in sequence, maintaining state"""
        channel = self.client.invoke_shell()
        output = ""
        
        for cmd in commands:
            channel.send(cmd + "\n")
            time.sleep(1.0)  # Give the command time to execute
            while channel.recv_ready():
                output += channel.recv(4096).decode()
        
        channel.close()
        return output

    def execute_command(self, command):
        """Execute a single command"""
        return self.execute_commands([command])

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.client:
            self.client.close()

# Update the connection functions to return the context manager
def connect_to_mmm():
    return SSHConnection(
        hostname="young.rc.ucl.ac.uk",
        username="mmm1486",
        key_path="/home/jakub/.ssh/id_rsa"
    )

def connect_to_hx1():
    return SSHConnection(
        hostname="login.hx1.hpc.ic.ac.uk",
        username="jl24018",
        password="1KhBRMU!VVgK7gzC4zHq"
    )


def execute_remote_command(ssh_client, command):
    """Execute command on remote system and return output"""
    stdin, stdout, stderr = ssh_client.exec_command(command)
    return stdout.read().decode(), stderr.read().decode()

def read_file(filepath):
    """Read CSV file using pandas"""
    df = pd.read_csv(filepath, sep=',')
    return df

def get_job_status(ssh_client, job_id):
    """Get simple status of a specific job"""
    command = f"qstat -j {job_id} 2>&1"  # 2>&1 captures both stdout and stderr
    output = ssh_client.execute_commands([command])
    
    if "Following jobs do not exist" in output:
        return "COMPLETED"
    elif "error" in output.lower():
        return "ERROR"
    else:
        # Get basic status using qstat
        output = ssh_client.execute_commands([f"qstat | grep {job_id}"])
        if output:
            print(output.split())
            # qstat output format: job-ID prior name user state
            return output.split()[3]  # Will return 'r' for running, 'qw' for queued, etc.
        return "UNKNOWN"


def submit_job_mmm(ssh_client, row):
    """Submit job to MMM"""
    # get info of job_id
    job_info = get_job_status(ssh_client, row['job_id'])

    print(job_info)
    assert False
    
    
    
    output = ssh_client.execute_commands([
        f"cd /home/mmm1486/projects/{project}/scripts/{experiment}",
        f"pwd"
    ])
    # qsub -N "$system" -v "FILEPATH=$filepath,SYSTEM=$system,OUTPUT_DIR=$output_dir,RESTART=$restart" $TEMPLATE
    print(output)

def submit_job_hx1(ssh_client, project, experiment):
    """Submit job to HX1"""
    pass

def submit_job(ssh_client, row):
    """Submit job to cluster"""
    if row['cluster'] == "mmm":
        submit_job_mmm(ssh_client, row)
    elif row['cluster'] == "hx1":
        submit_job_hx1(ssh_client, row)
    else:
        raise ValueError(f"Unknown cluster: {row['cluster']}")

if __name__ == "__main__":
    lines = read_file("jobs_to_submit.csv")
    
    for index, row in lines.iterrows():

        print(row['cluster'], row['project'], row['experiment'])
        # Use context manager for SSH connection
        with (connect_to_mmm() if row['cluster'] == "mmm" else connect_to_hx1()) as ssh_client:
            if ssh_client:
                submit_job(ssh_client, row)
        
        



        # jobs = get_jobs_status(ssh_client)
        # print(jobs)




