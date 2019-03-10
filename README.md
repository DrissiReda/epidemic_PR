# Epidemic Spreading Page Rank Simulation

A very simple epidemic simulation using pageRank algorithm. 

_Developed using Python **3.7.1**_

## Build

### Linux/ MacOS
- Download virtualenv
```sh
pip install virtualenv
```
- Create a virtualenv

``` sh
virtualenv venv
```

- Activate virtualenv

``` sh
virtualenv venv && source venv/bin/activate
```

- Install dependencies

``` sh
pip install -r requirements.txt
```
### Windows

- Install python3.5+ with tkinter support From [here](https://www.python.org/downloads/windows/)
- Install virtualenv

```PowerShell
pip install virtualenv virtualenvwrapper-powershell
mkdir '~\.virtualenvs'
```
- Activate virtualenv

```PowerShell
New-VirtualEnvironment venv
```

- Install dependencies

```PowerShell
pip install -r requirements.txt
```
## Run

``` sh
cd src
python -m app.main -g FILENAME
```
For more details you can check help with 
```sh
python -m app.main -h
```
