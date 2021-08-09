// beadplot.spec.js is created with Cypress, https://www.cypress.io/
// direct download of Cypress is required to run the tests
import clusters from '../fixtures/test_clusters.json'

describe('Array-input functions', () => {
    let coldates = []
    it('Unique', () => {
        cy.window().then((win) => {
            Object.values(clusters[0]["nodes"]).forEach(e => e.forEach(e => coldates.push(e[0])))
            expect(win.unique(coldates).length).to.eq(12)
            expect(win.unique(coldates)).to.eql(['2020-05-02','2020-05-15', '2020-05-26', '2020-06-10', '2020-06-16', 
            '2020-06-05', '2020-07-05', '2020-07-07', '2020-06-03', '2020-06-24', '2020-07-01', '2020-07-08'])
        })
    })
    it('Mode', () => {
        cy.window().then((win) => {
            expect(win.mode(coldates)).to.eq("2020-05-15")
        })
    })
    it('Tabulate', () => {
        cy.window().then((win) => {
            expect(win.tabulate(coldates)).to.eql({'2020-05-02': 1,'2020-05-15': 2, '2020-05-26': 1, '2020-06-10': 1, '2020-06-16': 1, 
            '2020-06-05': 1, '2020-07-05': 1, '2020-07-07': 1, '2020-06-03': 1, '2020-06-24': 1, '2020-07-01': 1, '2020-07-08': 1})
        })
    })
})

describe('Parse Functions', () => {
    var result = [];
    var beaddata; 
    it('parse_variant', () => {
        cy.window().then((win) => {
            let variant = Object.values(clusters[0].nodes)[0];
            let accn = Object.keys(clusters[0].nodes)[0];
            result = win.parse_variant(variant, 0, 0, accn, null, null); 

            expect(result.variant).to.eql({
                'accession': "EPI_ISL_542678",
                'label': "USA/TX-HMH-1371",
                'x1': new Date('2020-05-02 '),  
                'x2': new Date('2020-06-10 '),  
                'y1': 0,
                'y2': 0,
                'count': 5,
                'country': {"USA": 5},
                'region': ['North America','North America','North America','North America','North America'],
                'numBeads': 4,
                'parent': null,
                'dist': 0,
                'unsampled': false
            })
            expect(result.points[1]).to.eql({
                'cidx': 0,
                'variant': "EPI_ISL_542678",
                'x': new Date('2020-05-15 '),
                'y': 0,
                'count': 2,
                'accessions': ["EPI_ISL_545194", "EPI_ISL_545201"],
                'labels': ['USA/TX-HMH-2234', 'USA/TX-HMH-2246'],
                'region1': 'North America',
                'region': ['North America','North America'],
                'country': {"USA": 2},
                'parent': null,
                'dist': 0})
            })
    })
    it('parse_edgelist', () => {
        cy.window().then((win) => {
            // One variant 
            expect(win.parse_edgelist(clusters[0], [result.variant], [result.points])).to.eql([])

            // More than one variant
            let variants = [result.variant]
            variants.push({
                'accession': "EPI_ISL_543301",
                'label': "TX-HMH-4019",
                'x1': new Date('2020-06-16 '),  
                'x2': new Date('2020-06-16 '),  
                'y1': 0,
                'y2': 0,
                'count': 1,
                'country': {"USA": 1},
                'region': ['North America'],
                'numBeads': 1,
                'parent': null,
                'dist': 0,
                'unsampled': false
            })

            let points = [result.points]
            points = points.concat({
                'cidx': 1,
                'variant': "EPI_ISL_543301",
                'x': new Date("2020-06-16 "),
                'y': 0,
                'count': 1,
                'accessions': ["EPI_ISL_545450"],
                'labels': ['USA/TX-HMH-2649'],
                'region1': 'North America',
                'region': ['North America'],
                'country': {"USA": 1},
                'parent': null,
                'dist': 0
            })

            expect(win.parse_edgelist(clusters[0], variants, points)).to.eql([
                {
                    'y1': 0,
                    'y2': 0,
                    'x1': new Date("2020-06-16 "),  
                    'x2': new Date("2020-06-16 "),
                    'parent': "USA/TX-HMH-1371",
                    'child': "TX-HMH-4019",
                    'dist': 3.12,
                    'support': undefined
                }
            ])
        })
    })
    it('parse_clusters', () => {
        // Passing of the test is dependent on parse_edgelist and parse_variant functions
        cy.window().then((win) => {
            let edgelist = [{'y1': 1, 'y2': 2, 'x1': new Date("2020-06-16 "), 'x2': new Date("2020-06-16 "), 'parent': "USA/TX-HMH-1371", 'child': "USA/TX-HMH-4019", 'dist': 3.12, 'support': undefined},
            {'y1': 1, 'y2': 3, 'x1': new Date("2020-05-02 "), 'x2': new Date("2020-05-02 "), 'parent': "USA/TX-HMH-1371", 'child': "unsampled0", 'dist': 1.61, 'support': 0.62},
            {'y1': 1, 'y2': 4, 'x1': new Date("2020-06-05 "), 'x2': new Date("2020-06-05 "), 'parent': "USA/TX-HMH-1371", 'child': "USA/TX-HMH-2649", 'dist': 2.67, 'support': undefined},
            {'y1': 1, 'y2': 5, 'x1': new Date("2020-07-05 "), 'x2': new Date("2020-07-05 "), 'parent': "USA/TX-HMH-1371", 'child': "USA/TX-HMH-MCoV-10847", 'dist': 3.36, 'support': 0.94},
            {'y1': 3, 'y2': 6, 'x1': new Date("2020-06-03 "), 'x2': new Date("2020-06-03 "), 'parent': "unsampled0", 'child': "USA/TX-HMH-2767", 'dist': 1.24, 'support': undefined},
            {'y1': 3, 'y2': 7, 'x1': new Date("2020-06-24 "), 'x2': new Date("2020-06-24 "), 'parent': "unsampled0", 'child': "USA/TX-DSHS-1864", 'dist': 2.72, 'support': undefined},
            {'y1': 5, 'y2': 8, 'x1': new Date("2020-07-01 "), 'x2': new Date("2020-07-01 "), 'parent': "USA/TX-HMH-MCoV-10847", 'child': "USA/TX-HMH-MCoV-8639", 'dist': 0.87, 'support': undefined},
            {'y1': 5, 'y2': 9, 'x1': new Date("2020-07-08 "), 'x2': new Date("2020-07-08 "), 'parent': "USA/TX-HMH-MCoV-10847", 'child': "USA/TX-HMH-MCoV-10410", 'dist': 1.12, 'support': undefined}]

            let variants = [{"accession": "EPI_ISL_542678","label": "USA/TX-HMH-1371","x1": new Date("2020-05-02 "),"x2": new Date("2020-07-05 "),"y1": 1,"y2": 1,"count": 5,"country": {"USA": 5},
                    "region": [ "North America", "North America", "North America", "North America", "North America"],"numBeads": 4,"parent": null,"dist": 0,"unsampled": false },
                {"accession": "EPI_ISL_543301","label": "USA/TX-HMH-4019","x1": new Date("2020-06-16 "),"x2": new Date("2020-06-16 "),"y1": 2,"y2": 2,"count": 1,"country": {"USA": 1},
                    "region": [ "North America"],"numBeads": 1,"parent": "USA/TX-HMH-1371","dist": 3.12,"unsampled": false},
                {"accession": "unsampled0","label": "unsampled0","x1": new Date("2020-05-02 "),"x2": new Date("2020-07-08 "),"y1": 3,"y2": 3,"count": 0,"country": null,
                    "region": null,"numBeads": 0,"parent": "USA/TX-HMH-1371","dist": 1.61,"unsampled": true},
                {"accession": "EPI_ISL_545450","label": "USA/TX-HMH-2649","x1": new Date("2020-06-05 "),"x2": new Date("2020-06-05 "),"y1": 4,"y2": 4,"count": 1,"country": {"USA": 1},
                    "region": ["North America"],"numBeads": 1,"parent": "USA/TX-HMH-1371","dist": 2.67,"unsampled": false },
                {"accession": "EPI_ISL_789657","label": "USA/TX-HMH-MCoV-10847","x1": new Date("2020-07-01 "),"x2": new Date("2020-07-08 "),"y1": 5,"y2": 5,"count": 2,"country": {"USA": 2},
                    "region": ["North America","North America"],"numBeads": 2,"parent": "USA/TX-HMH-1371","dist": 3.36,"unsampled": false},
                {"accession": "EPI_ISL_545540","label": "USA/TX-HMH-2767","x1": new Date("2020-06-03 "),"x2": new Date("2020-06-03 "),"y1": 6,"y2": 6,"count": 1,"country": {"USA": 1},
                    "region": ["North America"],"numBeads": 1,"parent": "unsampled0","dist": 1.24,"unsampled": false},
                {"accession": "EPI_ISL_745299","label": "USA/TX-DSHS-1864","x1": new Date("2020-06-24 "),"x2": new Date("2020-06-24 "),"y1": 7,"y2": 7,"count": 1,"country": {    "USA": 1},
                    "region": ["North America"],"numBeads": 1,"parent": "unsampled0","dist": 2.72,"unsampled": false},
                {"accession": "EPI_ISL_780668","label": "USA/TX-HMH-MCoV-8639","x1": new Date("2020-07-01 "),"x2": new Date("2020-07-01 "),"y1": 8,"y2": 8,"count": 1,"country": {    "USA": 1},
                    "region": ["North America"],"numBeads": 1,"parent": "USA/TX-HMH-MCoV-10847","dist": 0.87,"unsampled": false},
                {"accession": "EPI_ISL_789458","label": "USA/TX-HMH-MCoV-10410","x1": new Date("2020-07-08 "),"x2": new Date("2020-07-08 "),"y1": 9,"y2": 9,"count": 1,"country": {"USA": 1},
                    "region": ["North America"],"numBeads": 1,"parent": "USA/TX-HMH-MCoV-10847","dist": 1.12,"unsampled": false}
            ]

            let points = [{"cidx": "0","variant": "EPI_ISL_542678","x": new Date("2020-05-02 "),"y": 1,"count": 1,"accessions": [
                        "EPI_ISL_542678"],"labels": ["USA/TX-HMH-1371"],"region1": "North America","region": ["North America"],"country": {"USA": 1},
                        "parent": null,"dist": 0},
                        {"cidx": "0","variant": "EPI_ISL_542678","x": new Date("2020-05-15 "),"y": 1,"count": 2,"accessions": 
                        ["EPI_ISL_545194","EPI_ISL_545201"],"labels": ["USA/TX-HMH-2234","USA/TX-HMH-2246"],"region1": "North America","region": ["North America","North America"],
                        "country": {"USA": 2},"parent": null,"dist": 0},
                        {"cidx": "0","variant": "EPI_ISL_542678","x": new Date("2020-05-26 "),"y": 1,"count": 1,"accessions": 
                        ["EPI_ISL_545290"],"labels": ["USA/TX-HMH-2403"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": null,"dist": 0},
                        {"cidx": "0","variant": "EPI_ISL_542678","x": new Date("2020-06-10 "),"y": 1,"count": 1,"accessions": 
                        ["EPI_ISL_543906"],"labels": ["USA/TX-HMH-2976"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1}, "parent": null,"dist": 0},
                        {"cidx": "0","variant": "EPI_ISL_543301","x": new Date("2020-06-16 "),"y": 2,"count": 1,"accessions": 
                        ["EPI_ISL_543301"],"labels": ["USA/TX-HMH-4019"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-1371","dist": 3.12},
                        {"cidx": "0","variant": "EPI_ISL_545450","x": new Date("2020-06-05 "),"y": 4,"count": 1,"accessions": 
                        ["EPI_ISL_545450"],"labels": ["USA/TX-HMH-2649"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-1371","dist": 2.67},
                        {"cidx": "0","variant": "EPI_ISL_789657","x": new Date("2020-07-05 "),"y": 5,"count": 1,"accessions": 
                        ["EPI_ISL_789657"],"labels": ["USA/TX-HMH-MCoV-10847"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-1371","dist": 3.36},
                        {"cidx": "0","variant": "EPI_ISL_789657","x": new Date("2020-07-07 "),"y": 5,"count": 1,"accessions": 
                        ["EPI_ISL_784276"],"labels": ["USA/TX-HMH-MCoV-8225"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-1371","dist": 3.36},
                        {"cidx": "0","variant": "EPI_ISL_545540","x": new Date("2020-06-03 "),"y": 6,"count": 1,"accessions": 
                        ["EPI_ISL_545540"],"labels": ["USA/TX-HMH-2767"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "unsampled0","dist": 1.24},
                        {"cidx": "0","variant": "EPI_ISL_745299","x": new Date("2020-06-24 "),"y": 7,"count": 1,"accessions": 
                        ["EPI_ISL_745299"],"labels": ["USA/TX-DSHS-1864"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "unsampled0","dist": 2.72},
                        {"cidx": "0","variant": "EPI_ISL_780668","x": new Date("2020-07-01 "),"y": 8,"count": 1,"accessions": 
                        ["EPI_ISL_780668"],"labels": ["USA/TX-HMH-MCoV-8639"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-MCoV-10847","dist": 0.87},
                        {"cidx": "0","variant": "EPI_ISL_789458","x": new Date("2020-07-08 "),"y": 9,"count": 1,"accessions": 
                        ["EPI_ISL_789458"],"labels": ["USA/TX-HMH-MCoV-10410"],"region1": "North America","region": ["North America"],
                        "country": {"USA": 1},"parent": "USA/TX-HMH-MCoV-10847","dist": 1.12}]

            beaddata = win.parse_clusters(clusters)
            expect(beaddata[0].edgelist).to.eql(edgelist)
            expect(beaddata[0].points).to.eql(points)
            expect(beaddata[0].variants).to.eql(variants)
                
        })
    })
    it('Merge tables ', () => {
        cy.window().then((win) => {
            expect(win.merge_tables(beaddata[0].variants.map(x => x.country))).to.eql({'USA': 13})
        })
    })
    it('Beadplot components', () => {
        cy.window().then((win) => {
            // Checks beadplot components
            var lg_beaddata = win.beaddata
            win.beaddata = beaddata

            win.beadplot(0)
            cy.get('.beadplot-content>svg>g').find('#USA-TX-HMH-1371').should('have.attr', 'stroke', '#777')
            cy.get('.beadplot-content>svg>g').find('#USA-TX-HMH-MCoV-10847').should('have.attr', 'stroke', '#777')
            cy.get('.beadplot-content>svg>g').find('#unsampled0').should('have.attr', 'stroke', '#ccc')
            cy.get('[stroke="#bbd"]').should('have.length', 4)

            cy.get('.beadplot-content>svg>g>text').should(($text) => {
                expect($text).to.have.length(9)
                expect($text.eq(0)).to.contain('USA/TX-HMH-1371')
                expect($text.eq(1)).to.contain('USA/TX-HMH-4019')
                expect($text.eq(2)).to.contain('unsampled0')
                expect($text.eq(3)).to.contain('USA/TX-HMH-2649')
                expect($text.eq(4)).to.contain('USA/TX-HMH-MCoV-10847')
                expect($text.eq(5)).to.contain('USA/TX-HMH-2767')
                expect($text.eq(6)).to.contain('USA/TX-DSHS-1864')
                expect($text.eq(7)).to.contain('USA/TX-HMH-MCoV-8639')
                expect($text.eq(8)).to.contain('USA/TX-HMH-MCoV-10410')
            }) 

            cy.get('.beadplot-content>svg>g>circle').then(($circ) => {
                expect($circ).to.have.length(12)
                cy.get($circ.eq(1)).should('have.attr', 'r', `${4*Math.sqrt(2)}`)
                for (let i = 0; i < 4; i++) {
                    cy.get($circ.eq(i)).should('have.attr', 'cy', 20)
                }

                cy.get($circ.eq(4)).should('have.attr', 'cy', 33.75)
                cy.get($circ.eq(5)).should('have.attr', 'cy', 61.25)
                cy.get($circ.eq(6)).should('have.attr', 'cy', 75)
                cy.get($circ.eq(7)).should('have.attr', 'cy', 75)
                cy.get($circ.eq(8)).should('have.attr', 'cy', 88.75)
                cy.get($circ.eq(9)).should('have.attr', 'cy', 102.5)
                cy.get($circ.eq(10)).should('have.attr', 'cy', 116.25)
                cy.get($circ.eq(11)).should('have.attr', 'cy', 130)
            })                 
        })
    })
})

describe('Edge slider', () => {
    var targetValue, currentValue, steps, arrows = 0
    const increment = 0.01

    it('Arrows increment slider by 0.01', () => {
        cy.get('#left-arrow').click()
        cy.get('#custom-handle').contains('1.99')

        cy.get('#right-arrow').click()
        cy.get('#custom-handle').contains('2')
    })
    it('Max slider value is equal to max edge length in cluster', () => {

        // https://stackoverflow.com/questions/64855669/moving-slider-with-cypress

        targetValue = 3.36
        currentValue = 2
        steps = (targetValue - currentValue) / increment + 1
        arrows = '{rightarrow}'.repeat(steps)

        cy.get('#custom-handle').contains(2).type(arrows)
        cy.get('#custom-handle').should('contain', 3.36)

        // Number displayed should not increase once max is reached
        cy.get('#custom-handle').type('{rightarrow}')
        cy.get('#custom-handle').should('contain', 3.36)
    })
    it('All edges are visible on max slider value ', () => {
        cy.get('[stroke="#bbd"]').should('have.length', 8)
    })
    it('3 edges are visible on slider value of 1.50 ', () => {
        targetValue = 1.50
        currentValue  = 3.36
        steps = (currentValue - targetValue) / increment + 1
        arrows = '{leftarrow}'.repeat(steps)

        cy.get('#custom-handle').type(arrows)
        cy.get('#custom-handle').should('contain', 1.50)
        cy.get('[stroke="#bbd"]').should('have.length', 3)
    })
    it('No edges are visible on slider value of 0.86 ', () => {
        targetValue = 0.86
        currentValue  = 1.50
        steps = (currentValue - targetValue) / increment
        arrows = '{leftarrow}'.repeat(steps)

        cy.get('#custom-handle').type(arrows)
        cy.get('#custom-handle').should('contain', 0.86)
        cy.get('[stroke="#bbd"]').should('not.exist')
    })
})

describe('Tooltips', () => {
    it('Parent, child, vertical edge and bead', () => {
        cy.get('[stroke-width="3"]').first().trigger('mouseover')
        cy.get('.tooltip').contains('North America: 5')
        cy.get('.tooltip').contains('Unique collection dates: 4')
        cy.get('.tooltip').contains('2020-05-02 / 2020-07-05')
 
        cy.get('[stroke="#ccc"]').first().trigger('mouseover')
        cy.get('.tooltip').contains('Parent: USA/TX-HMH-1371')
        cy.get('.tooltip').contains('Genomic distance: 1.61')

        cy.get('#right-arrow').click()
        cy.get('[stroke="#bbd"]').trigger('mouseover', {force: true})
        cy.get('.tooltip').contains('Parent: USA/TX-HMH-MCoV-10847')
        cy.get('.tooltip').contains('Child: USA/TX-HMH-MCoV-8639')
        cy.get('.tooltip').contains('Genomic distance: 0.87')
        cy.get('.tooltip').contains('Support: n/a')
        cy.get('.tooltip').contains('Collection date: 2020-07-01')

        cy.get('circle:visible').first().trigger('mouseover')
        cy.get('.tooltip').contains('North America: 1')
        cy.get('.tooltip').contains('Collection date: 2020-05-02')
    })
})

describe('Tables', () => {
    it('region_to_string', () => {
        cy.window().then((win) => {
            expect(win.region_to_string({'North America':5,'Asia':10})).to.eq("<b>Number of cases:</b><br>&nbsp;&nbsp;North America: 5<br>&nbsp;&nbsp;Asia: 10<br>Total: 15<br>")
        })
    })
    it('Country details: gentable', () => {
        cy.get('#tabs').contains('Countries').click()
        var tips = {'allregions': {'North America': 13}, 'country': {'Canada': 8, 'USA': 5}}

        cy.window().then((win) => {
            win.gentable(tips)
        })
        cy.get('table:visible>tbody>tr').should('have.length', 2)
        cy.get('table:visible>tbody>tr').eq(0).contains('North America')
        cy.get('table:visible>tbody>tr').eq(0).contains('Canada')
        cy.get('table:visible>tbody>tr').eq(0).contains('8')

        cy.get('table:visible>tbody>tr').eq(1).contains('North America')
        cy.get('table:visible>tbody>tr').eq(1).contains('USA')
        cy.get('table:visible>tbody>tr').eq(1).contains('5')
    })
    it('Sequence details: gen_details_table', () => {
        cy.get('#tabs').contains('Samples').click()

        cy.window().then((win) => {
            win.gen_details_table(win.beaddata[0].points[0])

            cy.get('#tabs-2').contains(`${win.beaddata[0].points[0].variant}`)
            cy.get('#tabs-2').contains(`${win.beaddata[0].points[0].labels[0]}`)
            cy.get('#tabs-2').contains(`${new Date(win.beaddata[0].points[0].x).toISOString().slice(0, 10)}`)
        })
    }) 
})