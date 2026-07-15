# ContextGuide

Classes that will use Encyclopedia's framework to implement Context with MProphet. This repo will serve as a source code for the ideas that are encompassed in Context.

## Installing ContextGuide from the command line

### Requirements

Before installing ContextGuide, confirm that the following programs are installed:

* Java Development Kit 17
* Apache Maven
* Git

Check the installations with:

```bash
java -version
javac -version
mvn -version
git --version
```

Java and Maven should report that they are using Java 17.

### Obtain the required dependency files

ContextGuide currently requires two dependencies that must be installed manually:

* Encyclopedia 6.6.24
* MSRawJava Core 26.7.1

You should receive the following four files from the ContextGuide maintainer:

```text
encyclopedia-6.6.24.jar
encyclopedia-6.6.24.pom
msrawjava-core-26.7.1.jar
msrawjava-core-26.7.1.pom
```

Place all four files in the same directory. For example:

```text
contextguide-dependencies/
├── encyclopedia-6.6.24.jar
├── encyclopedia-6.6.24.pom
├── msrawjava-core-26.7.1.jar
└── msrawjava-core-26.7.1.pom
```

Open a terminal and navigate to that directory:

```bash
cd /path/to/contextguide-dependencies
```

Replace `/path/to/contextguide-dependencies` with the actual directory containing the files.

### Install MSRawJava Core

Run:

```bash
mvn org.apache.maven.plugins:maven-install-plugin:3.1.4:install-file -Dfile=msrawjava-core-26.7.1.jar -DpomFile=msrawjava-core-26.7.1.pom
```

The command should end with:

```text
BUILD SUCCESS
```

### Install Encyclopedia

Run:

```bash
mvn org.apache.maven.plugins:maven-install-plugin:3.1.4:install-file -Dfile=encyclopedia-6.6.24.jar -DpomFile=encyclopedia-6.6.24.pom
```

The command should also end with:

```text
BUILD SUCCESS
```

These commands install the dependencies in the local Maven repository, normally located at:

```text
Linux, macOS, or WSL: ~/.m2/repository/
Windows:              %USERPROFILE%\.m2\repository\
```

The dependency JAR files should not be copied into the ContextGuide repository itself.

### Download ContextGuide

Navigate to the directory where you want to store the project, then clone the repository:

```bash
git clone https://github.com/shannon225/contextguide.git
```

Enter the repository:

```bash
cd contextguide
```

Use the `main` branch:

```bash
git switch main
git pull origin main
```

Confirm the active branch:

```bash
git branch --show-current
```

The command should print:

```text
main
```

### Compile ContextGuide

From the root of the ContextGuide repository, where `pom.xml` is located, run:

```bash
mvn clean compile
```

A successful compilation should end with:

```text
BUILD SUCCESS
```

Compiled class files will be placed under:

```text
target/classes/
```

### Verify the locally installed dependencies

To verify that Maven recognizes the dependencies, run:

```bash
mvn dependency:tree
```

The output should include:

```text
maccoss:encyclopedia:jar:6.6.24
org.searlelab:msrawjava-core:jar:26.7.1
```

### Updating ContextGuide later

To retrieve future updates:

```bash
cd contextguide
git switch main
git pull origin main
mvn clean compile
```

The two dependency files normally need to be installed only once for each computer and Maven installation.

## Troubleshooting

### Maven cannot find one of the dependency files

Confirm that the terminal is currently in the directory containing all four dependency files:

```bash
ls
```

On Windows PowerShell, use:

```powershell
Get-ChildItem
```

Check that the filenames in the installation command exactly match the actual filenames.

### Maven reports `Could not find artifact`

Confirm that the local repository contains the installed artifacts.

For Encyclopedia:

```bash
ls ~/.m2/repository/maccoss/encyclopedia/6.6.24/
```

For MSRawJava Core:

```bash
ls ~/.m2/repository/org/searlelab/msrawjava-core/26.7.1/
```

Each directory should contain both a `.jar` file and a `.pom` file.

### Maven created `.lastUpdated` files

A `.lastUpdated` file indicates that Maven previously failed to download or resolve an artifact.

After successfully installing the dependency locally, remove stale markers with:

```bash
find ~/.m2/repository/maccoss/encyclopedia/6.6.24 -name '*.lastUpdated' -delete
find ~/.m2/repository/org/searlelab/msrawjava-core/26.7.1 -name '*.lastUpdated' -delete
```

Then retry:

```bash
mvn clean compile
```

### Maven is using the wrong Java version

Run:

```bash
mvn -version
```

Confirm that the reported Java version is Java 17. If it is not, update `JAVA_HOME` so that Maven uses the Java 17 JDK.
