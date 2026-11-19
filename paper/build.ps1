param(
    [string]$Tectonic = 'tectonic',
    [string]$Bundle = 'https://data1b.fullyjustified.net/tlextras-2022.0r0.tar'
)
$ErrorActionPreference = 'Stop'
$build = Join-Path $PSScriptRoot 'build'
$public = Join-Path (Split-Path -Parent $PSScriptRoot) 'web/public/research'
[void][IO.Directory]::CreateDirectory($build)
[void][IO.Directory]::CreateDirectory($public)
Push-Location -LiteralPath $PSScriptRoot
try {
    & $Tectonic -b $Bundle --keep-logs --keep-intermediates --outdir $build manuscript.tex
    if ($LASTEXITCODE -ne 0) { throw "Tectonic exited with status $LASTEXITCODE." }
    Copy-Item -LiteralPath (Join-Path $build 'manuscript.pdf') -Destination (Join-Path $public 'biomass-supply-chain.pdf')
    Copy-Item -LiteralPath (Join-Path $PSScriptRoot 'manuscript.tex') -Destination $public
    Copy-Item -LiteralPath (Join-Path $PSScriptRoot 'references.bib') -Destination $public
} finally {
    Pop-Location
}
