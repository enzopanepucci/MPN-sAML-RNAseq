# Auditoria dos scripts recebidos

Os três scripts históricos foram preservados em `archive/`.

## O que foi confirmado

- A amostra `Amostra_18` corresponde ao segundo pré de P09 e foi excluída explicitamente.
- DESeq2 pareado: `design = ~ Patient_ID + Condition`, contraste `pos` versus `pre`.
- Filtro usado no script de origem: `rowSums(counts) > 30` e pelo menos cinco amostras com contagem não-zero.
- Ranking híbrido: `p < 0.05` usa `sign(log2FC) * (-log10(0.05) + abs(log2FC))`; `0.05 <= p < 0.1` usa `sign(log2FC) * (-log10(p) + abs(log2FC))`; `p >= 0.1` usa `sign(log2FC) * -log10(p)`.
- GSEA histórico usou Ensembl IDs, `fgsea`, KEGG com `minSize=15`, Hallmark com `minSize=10`, `maxSize=500` e seleção por `padj < 0.05`.

## O que não estava fechado

- A versão exata do R, `fgsea` e `msigdbr` usada na análise histórica.
- O arquivo original da matriz de assinatura do CIBERSORTx.
- Os outputs finais que deram exatamente os totais 98 KEGG e 47 Hallmark.

Por isso o projeto inclui a análise reproduzível, checkpoints e uma tabela de concordância com os números do manuscrito. A execução deve ser conferida antes de substituir qualquer número no artigo.
