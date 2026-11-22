param(
    [Parameter(Mandatory = $true)][string]$DocumentPath,
    [Parameter(Mandatory = $true)][string]$PdfPath
)

$ErrorActionPreference = 'Stop'
$documentFile = (Resolve-Path -LiteralPath $DocumentPath).Path
$pdfFile = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($PdfPath)
$word = New-Object -ComObject Word.Application
try {
    $word.Visible = $false
    $word.DisplayAlerts = 0
    $document = $word.Documents.Open($documentFile, $false, $true)
    try {
        # Export the complete document with heading bookmarks and structure tags.
        $document.ExportAsFixedFormat(
            $pdfFile, 17, $false, 0, 0, 1, 1, 0, $true, $false, 1, $true, $true, $false
        )
    } finally {
        $document.Close(0)
    }
} finally {
    $word.Quit(0)
}
