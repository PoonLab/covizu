// covizu.spec.js is created with Cypress, https://www.cypress.io/
// direct download of Cypress is required to run the tests

describe('CoVizu page', () => {
    it('Language options redirect', () => {
        cy.contains('a', 'es').click()
        cy.url().should('include', '/index-es.html')

        cy.contains('a', 'fr').click()
        cy.url().should('include', '/index-fr.html')

        cy.contains('a', 'zh').click()
        cy.url().should('include', '/index-zh.html')

        cy.contains('a', 'en').click()
        cy.url().should('include', '/index.html')
    })
    it('Welcome dialog opens automatically', () => {
        cy.get('#splash').should('be.visible')
        cy.get('#splash-button').should('be.disabled')
    })
    it('Splash button is enabled after clusters.json loads', () => {
        cy.get('#splash-button', {timeout: 40000}).should('be.enabled')
        cy.get('#splash-button').click().should('not.be.visible')
    })
})

describe('Help dialogs', () => {
    it('Open and close', () => {
        cy.get('[style="cursor: help"]').each(($el) => {
            cy.get($el).click()
            cy.get('.ui-dialog:visible').find('.ui-dialog-titlebar-close').click()
            cy.get('.ui-dialog').should('not.be.visible')
        })
    })
})

describe('Search Interface', () => {
    it('Buttons are initially disabled', () => {
        cy.get('#navigation>button').each(($el) => {
            cy.get($el).should('be.disabled')
        })
    })
    it('Error message for invalid date formats', () => {
        cy.get('#start-date').type('2021-05-13')
        cy.get('#end-date').type('2021-05-10')
        cy.get('#search-button').click({force:true})
        cy.contains('Start Date must be before the End Date.')

        cy.get('#start-date').clear()
        cy.get('#end-date').clear()
        cy.get('#start-date').type('20210510')
        cy.get('#search-button').click({force:true})
        cy.contains('Invalid Start Date (YYYY-MM-DD)')

        cy.get('#start-date').clear()
        cy.get('#end-date').clear()
        cy.get('#end-date').type('20210513')
        cy.get('#search-button').click({force:true})
        cy.contains('Invalid End Date (YYYY-MM-DD)')
        cy.get('#end-date').clear()
    })
    it("Searching 'HMH' results in selection and expected point count", () => {
        cy.get('#search-input').type('HMH')
        cy.get('#search-button').click({force:true})
        cy.get('.selectionH').should('exist')
        cy.get('.SelectedCluster').should('exist')
        cy.window().then((win) => {
            cy.get('#search_stats').contains(`1 of ${win.bead_id_to_accession.length} points`)
        })
    })
    it("Clicking the next button highlights the correct bead", () => {
        cy.window().then((win) => {
            cy.get('.selectionH').should('have.attr', 'bead', `${win.bead_id_to_accession[0]}`)
            cy.get('#next_button').click()
            cy.get('.selectionH').should('have.length', 1)
            cy.get('.selectionH').should('have.attr', 'bead', `${win.bead_id_to_accession[1]}`)
        })
    })
    it('Clear button resets search interface and cluster selection', () => {
        cy.get('#clear_button').click()
        cy.get('#search-bar>input').each(($el) => {
            cy.get($el).should('be.empty')
        })
        cy.get('#search_stats').contains('0 of 0 points')
        cy.get('.selectionH').should('not.exist')
        cy.get('.SelectedCluster').should('not.exist')
    })
})

describe('Tooltips', () => {
    it('Appear on hover over cluster', () => {
        cy.get('rect:visible').first().trigger('mouseover')
        cy.get('.tooltip').should('be.visible').should('have.length', 1)
    })
    it('Cluster tooltips contain relevant and correct information', () => {
        cy.get('rect:visible').first().trigger('mouseover').invoke('attr','id').as('id')
        cy.get('@id').then(id => {
            var id_number = id.substring(3)
            cy.window().then((win) => {
                // Number of variants
                cy.get('.tooltip').contains(`Number of variants: ${win.tips[id_number]['varcount']}`)
                // Mutations list
                cy.wrap(win.tips[id_number]['mutations']).each(($el) => {
                    cy.get('.tooltip').contains($el)
                })
                // Regions: Num of Cases 
                cy.wrap(Object.keys(win.tips[id_number]['allregions'])).each(($el, index) => {
                    cy.get('.tooltip').contains(`${$el}: ${Object.values(win.tips[id_number]['allregions'])[index]}`)
                })
                // Mean diffs from root
                cy.get('.tooltip').contains(`Mean diffs from root: ${win.tips[id_number]['mean_ndiffs'].toFixed(2)}`) 
                // Deviation from clock
                cy.get('.tooltip').contains(`Deviation from clock: ${win.tips[id_number]['residual'].toFixed(2)}`)
                // Collection dates (UTC)
                cy.get('.tooltip').contains(`${win.tips[id_number]['first_date'].toISOString().slice(0, 10)} / ${win.tips[id_number]['last_date'].toISOString().slice(0, 10)}`)
            })
        })
    })
    it('Appear on hover over bead', () => {
        cy.get('circle:visible').first().trigger('mouseover')
        cy.get('.tooltip').should('be.visible').should('have.length', 1)
    })
    it('Beadplot tooltips contain relevant information', () => {
        cy.get('rect:visible').first().click().invoke('attr','id').as('lineage_id')
        cy.get('@lineage_id').then(id => {
            var lineage_id = id.substring(3)
            cy.window().then((win) => {
                // Tooltip for first line with stroke-width of 3 (Assuming that it is parent)
                cy.get('[stroke-width="3"]').first().trigger('mouseover', {force: true})
                cy.get('.tooltip').contains(`Unique collection dates: ${win.beaddata[win.tips[lineage_id]["cluster_idx"]]["variants"][0]["numBeads"]}`)
                cy.get('.tooltip').contains(`${win.beaddata[win.tips[lineage_id]["cluster_idx"]]["variants"][0]["x1"].toISOString().slice(0, 10)} / ${win.beaddata[win.tips[lineage_id]["cluster_idx"]]["variants"][0]["x2"].toISOString().slice(0, 10)}`)
                var regions = win.tabulate(win.beaddata[win.tips[lineage_id]["cluster_idx"]]["variants"][0]["region"])
                cy.wrap(Object.keys(regions)).each(($el, index) => {
                    cy.get('.tooltip').contains(`${$el}: ${Object.values(regions)[index]}`)
                })
            })
        })
    })
})