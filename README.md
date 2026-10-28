# MammoStationarity

**Análise da estacionariedade do ruído em imagens médicas utilizando Python.**

## Sobre o projeto

O MammoStationarity é uma ferramenta computacional desenvolvida em Python para estudar a distribuição espacial do ruído em imagens médicas, com foco inicial em imagens de mamografia.

O programa realiza uma varredura da imagem por regiões de interesse (*Regions of Interest* — ROIs), calcula o desvio padrão dos valores dos pixels em cada região e gera mapas de calor para visualizar como a variabilidade do sinal se distribui espacialmente.

A proposta é transformar imagens médicas em dados quantitativos, facilitando a investigação de padrões e a avaliação da estacionariedade do ruído.

## Objetivos

- Investigar a distribuição espacial do ruído em imagens médicas.
- Automatizar o cálculo de métricas estatísticas em diferentes regiões da imagem.
- Visualizar a variabilidade espacial por meio de mapas de calor.
- Apoiar estudos de qualidade de imagem e análise de dados aplicados à Física Médica.

## Como funciona

O programa utiliza uma abordagem baseada em regiões de interesse móveis para analisar a imagem de maneira sistemática.

1. **Leitura da imagem:** a imagem é carregada para análise computacional.
2. **Varredura por ROIs:** regiões de interesse percorrem a imagem seguindo um intervalo definido entre posições consecutivas.
3. **Cálculo estatístico:** o desvio padrão dos valores dos pixels é calculado em cada região.
4. **Visualização:** os resultados são organizados espacialmente e representados por mapas de calor, permitindo observar variações locais e possíveis padrões na distribuição do ruído.

### Parâmetros da análise

- **Tamanho da ROI:** 5 × 5 mm².
- **Passo da varredura:** 10 pixels.
- **Métrica estatística:** desvio padrão dos valores dos pixels.
- **Representação dos resultados:** mapas de calor da variabilidade espacial.

Os parâmetros podem ser ajustados conforme as características da imagem, os objetivos da análise e o protocolo que será seguido.

## Tecnologias utilizadas

- **Python** — desenvolvimento e automação da análise.
- **Estatística** — cálculo e interpretação de métricas de variabilidade.
- **Processamento de imagens** — varredura sistemática e análise de regiões da imagem.
- **Visualização de dados** — representação espacial dos resultados por mapas de calor.

## Aplicações

O MammoStationarity foi desenvolvido para apoiar estudos de estacionariedade do ruído em imagens médicas, especialmente em mamografia digital.

A abordagem também pode ser investigada em outras modalidades de imagem, desde que sejam consideradas as características dos dados e as adaptações metodológicas necessárias.

## Contexto acadêmico

Este projeto está relacionado à pesquisa em Física Médica e ao estudo de métodos quantitativos para avaliação da qualidade de imagens médicas.

Seu desenvolvimento reúne programação, análise estatística e visualização de dados para investigar um problema experimental por meio de uma ferramenta computacional.

## Limitações

Os resultados dependem das características da imagem, dos parâmetros utilizados e das condições da análise. Os mapas de calor devem ser interpretados considerando o contexto físico e estatístico do problema.

O programa tem finalidade de pesquisa e análise experimental, não se destinando, por si só, ao diagnóstico clínico ou à tomada de decisões médicas.

## Próximos passos

- Aprimorar a documentação e a reprodutibilidade das análises.
- Expandir a avaliação para diferentes conjuntos de imagens.
- Investigar métodos complementares de análise estatística.
- Melhorar a apresentação e a comparação dos resultados.

---

**Autor:** Lucas Barbosa Araújo da Silva

**Área:** Física Médica · Python · Análise de Dados · Processamento de Imagens
