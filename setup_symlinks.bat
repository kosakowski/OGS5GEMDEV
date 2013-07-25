::: Symlink ExampleData :::
mkdir sources\Build\bin\Release\ExampleData
mkdir sources\Build\bin\Debug\ExampleData
Libs\ToolsWin\junction.exe .\sources\Build\bin\Release\ExampleData .\ExampleData
Libs\ToolsWin\junction.exe .\sources\Build\bin\Debug\ExampleData .\ExampleData