---
sidebar_position: 1
---




# Tutorial Intro

import useBaseUrl from '@docusaurus/useBaseUrl';


Let's discover **MontyCarlo in less than 5 minutes**. Just be sure you've installed it first! (The installation will take a lot more than 5 min...)

## Getting Started

Get started by **creating a new project**. Go to a terminal and activate the virtual environment that holds your MontyCarlo installation. To open the MontyCarlo shell, type the command

```shell
$ myco
```

This shell supports the usual navigation commands (ls and cd). Apart from that, the most important commands are:

- `create [name of your cool project]`: this creates a new project folder
- `open [name of your cool project]`: this opens the project
- `build`: this takes your source code and compiles it to pure Python
- `run`: runs the simulation !

The `create` command already writes some boilerplate code for you. 

## Running the first simulation

1. Open the MontyCarlo shell

```
$ myco
```

2. Create a project

```
myco@users/user_name/projects> create my_cool_project
```

3. Open your project

```
myco@users/user_name/projects> open my_cool_project
```

4. Build the project

```
myco@my_cool_project> build
```

5. Run the simulation

```
myco@my_cool_project> run
```

In the end it should look something like this:

<img alt="Docusaurus with Keytar" src={useBaseUrl('/img/tutorial/shell.png')} />


and the result of the simulation should be the following plot

<img alt="Docusaurus with Keytar" src={useBaseUrl('/img/tutorial/first_run.png')} />
